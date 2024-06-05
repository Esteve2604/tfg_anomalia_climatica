#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>
# include <thrust/sort.h>
# include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#define FILE_NAME_PREC "../copernicus_prec_data_clean_bien.nc"
#define FILE_NAME_TG "../copernicus_tg_data_clean_bien.nc"
#define NDIMS 3
#define LAT 434
#define LON 575
#define TIME 5844
#define NUM_THREADS 20
#define ERRCODE 2
#define MEMERRORCODE 3
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define MEMERR(e) {printf("Error: %s\n", e); exit(MEMERRORCODE);}

#define SHORT 1
#define DOUBLE 0
//Devuelve la longitud del dia
double getDayLength(double day);

//Devuelve el drought code calculado
double calculateDroughtCode(short temperature, double rainfall, double droughtCode, double day);

//Deja los datos calculados en estres_termico y ordena de menor a mayor part_tg_val
int calculateEstresTermico(short *part_tg_val, short *orderedTg, short *estres_termico);

//Deja los datos calculados en estres_hidrico y ordena de menor a mayor part_prec_val
int calculateEstresHidrico(double *part_prec_val, double *orderedPrec, short *estres_hidrico);

//Deja los datos calculados en anomalia_climatica
int calculateAnomaliaClimatica(short *estres_termico, short *estres_hidrico, short* anomalia_climatica);

//Función para encontrar el percentil al que pertenece un valor, devuelve el percentil
int encontrarPercentilShort(short *percentiles, short val);
int encontrarPercentilDouble(double *percentiles, double val);
int ordenarThrustGpu(short *);
//En esta versión se intenta calcular concurrentemente todos los LON de una LAT
int main(){
  int retval;
  int i,j,k, ncIdPrec[NUM_THREADS], ncIdTg[NUM_THREADS], threadId, total_pos;
  size_t *startp, *countp;
  short *part_tg_val, *device_part_tg_val, *estres_hidrico, *estres_termico, *anomalia_climatica, *orderedTg;
  double *part_prec_val, *drought_code, last_drought_code, *orderedPrec;
  int  despl_lat, despl_total;
  int precId, tgId, timeId;
  int invalidValue, notCalculated=0;
  double time[TIME];
  struct timeval t, t2;
  double segundos;
  gettimeofday(&t, NULL);
  total_pos = LAT*LON;
  printf("**Inicio** \n");
  printf("**Abriendo ficheros** \n");
  for(i=0; i<NUM_THREADS; i++){
    if((retval= nc_open(FILE_NAME_PREC, NC_NOWRITE, &ncIdPrec[i]))) //Abrimos fichero precipitacion
      ERR(retval);
    if((retval=nc_open(FILE_NAME_TG, NC_NOWRITE, &ncIdTg[i]))) //Abrimos fichero temperatura
      ERR(retval);
  }
  nc_get_var_double(ncIdPrec[0], timeId, time); //Inicializamos los dias
#pragma omp parallel private(retval,i,j,k,startp, countp, part_tg_val, device_part_tg_val, estres_hidrico, estres_termico,anomalia_climatica, part_prec_val, drought_code,last_drought_code,despl_lat,despl_total,invalidValue,timeId, orderedPrec,orderedTg, notCalculated, precId,tgId, threadId) shared(ncIdPrec,ncIdTg,time)  num_threads(NUM_THREADS)
  {
    threadId = omp_get_thread_num();
    if((retval=nc_inq_varid(ncIdPrec[threadId], "rr", &precId))) //Cogemos el id de la variable precipitacion
      ERR(retval);
    if((retval=nc_inq_varid(ncIdTg[threadId], "tg", &tgId))) //Cogemos el id de la variable temperatura
      ERR(retval);
    if((retval=nc_inq_varid(ncIdPrec[threadId], "time", &timeId))) //Cogemos el id de la variable tiempo del fichero de prec, aunque debería coincidir con el del fichero de tg.
      ERR(retval);
    printf("**Reservando memoria** \n");
    if((orderedPrec = (double *)malloc(TIME*sizeof(part_prec_val)*LON))==NULL)
      MEMERR("Reserva orderedPrec");
    if((orderedTg = (short *)malloc(TIME*sizeof(part_tg_val)*LON)) == NULL)
      MEMERR("Reserva orderedTg");
    if((part_tg_val=(short *)calloc(sizeof(short), TIME*LON)) ==NULL)
      MEMERR("Reserva part_tg_val");
    cudaMalloc((void**) &device_part_tg_val, TIME * sizeof(short)*LON);
    if((part_prec_val=(double *)calloc(sizeof(double),TIME*LON)) ==NULL)
      MEMERR("Reserva part_prec_val");
    if((drought_code = (double *)calloc(sizeof(double),TIME*LON))==NULL)
      MEMERR("Reserva drought_code");
    if((estres_termico = (short *)calloc(sizeof(short),TIME*LON))==NULL)
      MEMERR("Reserva estres_termico");
    if((estres_hidrico = (short *)calloc(sizeof(short),TIME*LON))==NULL)
      MEMERR("Reserva estres_hidrico");
    if((anomalia_climatica = (short *)calloc(sizeof(short),TIME*LON))==NULL)
      MEMERR("Reserva anomalia_climatica");
    startp=(size_t *)malloc(NDIMS*sizeof(size_t));
    countp=(size_t *)malloc(NDIMS*sizeof(size_t));
    *countp=1; //Leemos un dato de lat
    *(countp+1)=LON; //tantos datos como LON
    *(countp+2)=TIME; //y tantos datos como dias
    *(startp+2)=0;//Siempre empezamos a leer desde el primer dia
    i=0;
    printf("**Iniciando calculos proceso: %d** \n", omp_get_thread_num());
#pragma omp for schedule(guided)
    for(j=0;j<LAT; j++){
      if(j == LAT/4)
	printf("**25%% del calculo completado**\n");
      else if(j == LAT/2)
	printf("**50%% del calculo completado**\n");
      else if(j == ((LAT/4)*3))
	printf("**75%% del calculo completado**\n");
      despl_lat = j * LON; //Desplazamiento Filas
      invalidValue=0;
      despl_total = despl_lat + k; //Desplazamiento total con desplazamiento de columna
      //Obtención prec y tg del lugar de todos los días.
      *startp = j;
      *(startp+1)= 0;//Siempre empezamos en la primera LON de la lat
	
      if((retval=nc_get_vara_short(ncIdTg[threadId],tgId,startp,countp, part_tg_val))){
	printf("error get tg\n");
	ERR(retval);
      }
      cudaMemcpy(device_part_tg_val, part_tg_val,LON * TIME * sizeof(short),
		 cudaMemcpyHostToDevice);
      ordenarThrustGpu(device_part_tg_val);
      if((retval=nc_get_vara_double(ncIdPrec[threadId],precId,startp,countp, part_prec_val))){
	printf("error get prec\n");
	ERR(retval);
      }

      for(k=0;i<LON;k++){
	//Calculo del drought code
	last_drought_code=calculateDroughtCode(*(part_tg_val+k*LON), *(part_prec_val+k*LON), 0, time[0]);
      
	if(last_drought_code<0){
	  invalidValue = 1;
	} else {
	  drought_code[0]=last_drought_code;
	  for(i=1;i<TIME;i++){
	    last_drought_code=calculateDroughtCode(*(part_tg_val+k*LON+i), *(part_prec_val+k*LON+i), last_drought_code,time[i]);
	    *(drought_code+k*LON+i)=last_drought_code;
	    if(last_drought_code<0){
	      invalidValue = 1;
	      break;
	    }	
	  }
	}
      }
    
      //if(!invalidValue){
      //Cálculo del estrés hídrico de cada día
      //Ordenamos los datos para calcular los percentiles
      memcpy(orderedPrec, part_prec_val, TIME);
      for(i=0;i<LON;i++){
	thrust::sort(thrust::host, orderedPrec+i*LON, orderedPrec+i*LON + TIME);
	calculateEstresHidrico(part_prec_val+i*LON,orderedPrec+i*LON, estres_hidrico+i*LON);
      }
      //Ordenamos para el cálculo de los percentiles
      cudaMemcpy(orderedTg, device_part_tg_val, TIME * sizeof(short),
		 cudaMemcpyDeviceToHost);
      //Cálculo del estrés térmico de cada día
      for(i=0;i<LON;i++){
	calculateEstresTermico(part_tg_val+i*LON, orderedTg+i*LON, estres_termico+i*LON);
	//Cálculo de la anomalía climática
	calculateAnomaliaClimatica(estres_termico+i*LON, estres_hidrico+i*LON, anomalia_climatica+i*LON);
      }
      /* } else {
      //printf("Invalid values in LAT: %d LON: %d",j,k);
      notCalculated++;
      for(i=0; i<TIME; i++){
      *(estres_termico+i)=-1;
      *(estres_hidrico+i)=-1;
      *(anomalia_climatica+i) = -1;
      }
      }*/
    }
    printf("**Cálculos finalizados con éxito proceso: %d** \n", omp_get_thread_num());

    free(orderedTg);
    free(orderedPrec);
    free(startp);
    free(countp);
    free(drought_code);
    free(part_prec_val);
    free(part_tg_val);
    free(estres_termico);
    free(estres_hidrico);
    free(anomalia_climatica);
    cudaFree(device_part_tg_val);
    nc_close(ncIdPrec[threadId]);
    nc_close(ncIdTg[threadId]);
  }
 
  gettimeofday(&t2, NULL);
  segundos = (((t2.tv_usec - t.tv_usec)/1000000.0f)  + (t2.tv_sec - t.tv_sec));

  printf("\n *******> Duracion total: %f segundos\n\n", segundos);
  
  printf("Total:%lu , NotCalculated: \n",(unsigned long)LAT*LON*TIME);
  return 0;
}
int ordenarThrustGpu(short *device_part_tg_val){
 //Auxiliar que guarda a que posición pertenece cada elementos
  thrust::device_vector<float> d_keys(LON * TIME);

    // --- Generate row indices: Las genera creando indices desde 0 hasta SizeToOrder * NUM_ITER, después crea un iterador con el SizeToOrder constante, y divide, de forma que queda solo el NUM_ITER correspondiente a cada elemento, es decir, la posición a la que pertenece e inicializa en d_keys
    thrust::transform(thrust::make_counting_iterator(0),
                      thrust::make_counting_iterator(LON * TIME),
                      thrust::make_constant_iterator(LON),
                      d_keys.begin(),
                      thrust::divides<int>());

    // --- Back-to-back approach
    //Ordenamos todos los elementos del vector
    thrust::stable_sort_by_key(device_part_tg_val,
                               device_part_tg_val + LON * TIME,
                               d_keys.begin(),
                               thrust::less<short>());
    //Ordenamos las keys. De esta forma, al estar ordenado todo el vector, al ordenar las keys, cada posición estará ordenada
    // -6 2 4 1 3 5 ==> -6 1 2 3 4 5 ==> -6 2 4 1 3 5
    //  0 0 0 1 1 1 ==>  0 1 0 1 0 1 ==>  0 0 0 1 1 1
    thrust::stable_sort_by_key(d_keys.begin(),
                               d_keys.end(),
                               device_part_tg_val,
                               thrust::less<int>());
  return 0;
}
double calculateDroughtCode(short temperature, double rainfall, double droughtCode, double day){
  double prevDroughtCode= droughtCode;
  double efRainfall, moisture, prevMoisture, evapotranspiration;
  if(rainfall < 0 || temperature < -8000 || temperature > 8000 || droughtCode < 0){
    return -1;
  }
  if(rainfall>2.8){
    efRainfall = 0.86 * rainfall -1.27;
    prevMoisture = 800 * exp(-droughtCode/400);
    moisture = prevMoisture + 3.937 * efRainfall;
    prevDroughtCode = 400 * log(800/moisture);
    prevDroughtCode = prevDroughtCode < 0 ? 0 : prevDroughtCode;
  }
  evapotranspiration = 0.36 * (temperature + 2.0) + getDayLength(day);
  evapotranspiration = evapotranspiration < 0 ? 0 : evapotranspiration;
  return prevDroughtCode + 0.5 * evapotranspiration;
}

double getDayLength(double day){
  int dayYear = (int) day%365;
  if(dayYear < 90)
    return -1.6;
  if(dayYear<120)
    return 0.9;
  if(dayYear<151)
    return 3.8;
  if(dayYear<181)
    return 5.8;
  if(dayYear<212)
    return 6.4;
  if(dayYear<243)
    return 5.0;
  if(dayYear<273)
    return 2.4;
  if(dayYear<304)
    return 0.4;
  if(dayYear<365)
    return -1.6;
  return -9999;
}

int calculateEstresTermico(short *part_tg_val, short *orderedTg, short *estres_termico){
  int i;
  int tamano_percentil = TIME/100;
  short *percentiles = (short *) malloc(99*sizeof(short));
  short media_ocho_dias;
  //Sacamos los números que dividen los percentiles
  for(i=1; i<99; i++)
    percentiles[i-1]=*(orderedTg+i*tamano_percentil);
  for(i=0; i<7; i++) //Se necesita la media de los 8 días anteriores
    *(estres_termico) =-1;
  for(i=7; i<TIME; i++){
    media_ocho_dias = (short) (((int) *(part_tg_val+i-7) + (int) *(part_tg_val+i-6) + (int) *(part_tg_val+i-5) + (int) *(part_tg_val+i-4) + (int) *(part_tg_val+i-3) + (int) *(part_tg_val+i-2) + (int) *(part_tg_val+i-1) + (int) *(part_tg_val+i)) /8);
    *(estres_termico) = encontrarPercentilShort(percentiles,media_ocho_dias);
  }
  free(percentiles);
  return 0;
}

int calculateEstresHidrico(double *part_prec_val,double *orderedPrec, short *estres_hidrico){
  int i;
  int tamano_percentil = TIME/100;
  double *percentiles = (double *)malloc(99*sizeof(double));
  //Sacamos los números que dividen los percentiles
  for(i=1; i<99; i++)
    percentiles[i-1]=*(orderedPrec+i*tamano_percentil);
  //Sacamos los números que dividen los percentiles
  for(i=1; i<99; i++)
    percentiles[i-1]=*(orderedPrec+i*tamano_percentil);
  for(i=0; i<TIME; i++){ //Se necesita la media de los 8 días anteriores
    *(estres_hidrico) = encontrarPercentilDouble(percentiles, *(part_prec_val+i));
  }
  free(percentiles);
  return 0;
}

int calculateAnomaliaClimatica(short *estres_termico, short *estres_hidrico, short* anomalia_climatica){
  int i = 0;
  for(i=0; i<TIME; i++)
    *(anomalia_climatica+i) = (*(estres_termico+i) + *(estres_hidrico+i)) / 2;
  return 0;
}

int encontrarPercentilShort(short *percentiles, short val){
  int last_percentil=-1, first=0, middle=49, last=99;
  while(first<=last){
    if(*(percentiles+middle) < val){
      last_percentil=middle;
      first = middle +1;
    } else if(*(percentiles+middle)== val){
      last_percentil = middle;
      break;
    } else {
      last = middle -1;
    }
    middle = (first +last) /2;
  }
  return last_percentil+1; //+1 porque el percentil empieza en 0
}
int encontrarPercentilDouble(double *percentiles, double val){
  int last_percentil=-1, first=0, middle=49, last=99;
  while(first<=last){
    if(*(percentiles+middle) < val){
      last_percentil=middle;
      first = middle +1;
    } else if(*(percentiles+middle)== val){
      last_percentil = middle;
      break;
    } else {
      last = middle -1;
    }
    middle = (first +last) /2;
  }
  return last_percentil +1; //+1 porque el percentil empieza en 0
}
