#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>
#include <thrust/sort.h>
#include <thrust/async/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/async/copy.h>
#include <thrust/async/transform.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#define FILE_NAME_PREC "../copernicus_prec_data_clean.nc"
#define FILE_NAME_TG "../copernicus_tg_data_clean.nc"
#define NDIMS 3
#define LAT 437
#define LON 592
#define TIME 5844
#define NUM_THREADS 1
#define ERRCODE 2
#define MEMERRORCODE 3
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define MEMERR(e) {printf("Error: %s\n", e); exit(MEMERRORCODE);}

#define SHORT 1
#define DOUBLE 0
//Devuelve la longitud del dia
double getDayLength(double day);

//Devuelve el drought code calculado
int calculateDroughtCode(short *temperature, float *rainfall, float *droughtCode, double *days);

//Deja los datos calculados en estres_termico y ordena de menor a mayor part_tg_val
int calculateEstresTermico(short *part_tg_val, short *orderedTg, short *estres_termico);

//Deja los datos calculados en estres_hidrico y ordena de menor a mayor part_prec_val
int calculateEstresHidrico(float *drought_code, float *orderedDroughtCode, short *estres_hidrico);

//Deja los datos calculados en anomalia_climatica
int calculateAnomaliaClimatica(short *estres_termico, short *estres_hidrico, short* anomalia_climatica);

//Función para encontrar el percentil al que pertenece un valor, devuelve el percentil
short encontrarPercentilShort(short *percentiles, short val);
short encontrarPercentilFloat(float *percentiles, float val);
//En esta versión se intenta calcular concurrentemente todos los LON de una LAT
int main(){
  int retval;
  double *lon, *lat;
  int i,j,k, ncIdPrec[NUM_THREADS], ncIdTg[NUM_THREADS], threadId, total_pos;
  size_t *startp, *countp;
  int newNcId, newLonDim, newLatDim, newTimeDim, dimIds[3], newAnomaliaClimaticaId, newLonId, newLatId, newTimeId; //Ids for new file
  short *part_tg_val, *device_part_tg_val[NUM_THREADS], *estres_hidrico, *estres_termico, *anomalia_climatica, *orderedTg;
  float *part_prec_val,*device_drought_code, *drought_code, last_drought_code, *orderedDroughtCode;
  short *d_keys[NUM_THREADS]; //Para back to back
  int  despl_lat, despl_total, invalidValue[LON];
  int precId, tgId, timeId, lonId, latId;
  unsigned long notCalculated=0;
  size_t writeStartp[3], writeCountp[3];
  cudaStream_t streams[NUM_THREADS];
  cudaError_t err;
  double time[TIME];
  struct timeval t, t2;
  double segundos;
  gettimeofday(&t, NULL);
  total_pos = LAT*LON;
  printf("**Inicio** \n");
  if((lat = (double *)malloc(LAT*sizeof(double)))==NULL)
      MEMERR("Reserva lon");
  if((lon = (double *)malloc(LON*sizeof(double)))==NULL)
      MEMERR("Reserva lat");
  printf("**Abriendo ficheros** \n");
  for(i=0; i<NUM_THREADS; i++){
    if((retval= nc_open(FILE_NAME_PREC, NC_NOWRITE, &ncIdPrec[i]))) //Abrimos fichero precipitacion
      ERR(retval);
    if((retval=nc_open(FILE_NAME_TG, NC_NOWRITE, &ncIdTg[i]))) //Abrimos fichero temperatura
      ERR(retval);
    // Create the stream:
    err = cudaStreamCreate(&streams[i]);
    if (err != cudaSuccess)
      {
	printf("Error creating stream: %s\n", cudaGetErrorString(err));
	exit(1);
      }
    if((cudaMalloc((void**) &device_part_tg_val[i], TIME * sizeof(short)*LON)))
      MEMERR("Reserva device part tg val");
    if((cudaMalloc((void**) &d_keys[i], TIME * sizeof(short)*LON)))
      MEMERR("Reserva device part tg val");
  }
  printf("**Leyendo variables basicas** \n");
  if((retval=nc_inq_varid(ncIdPrec[0], "time", &timeId))) //Cogemos el id de la variable tiempo del fichero de prec, aunque debería coincidir con el del fichero de tg.
    ERR(retval);
  nc_get_var_double(ncIdPrec[0], timeId, time); //Inicializamos los dias
  if((retval=nc_inq_varid(ncIdPrec[0], "longitude", &lonId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdPrec[0], lonId, lon)))
    ERR(retval);
  if((retval=  nc_inq_varid(ncIdPrec[0], "latitude", &latId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdPrec[0], latId, lat)))
    ERR(retval);
  //Creamos el fichero final de anomalía climatica
   printf("**Creando fichero de anomalia climatica** \n");
  if((retval=nc_create("copernicus_anomalia_climatica_1995_2010.nc", NC_CLOBBER, &newNcId)))
    ERR(retval);
  printf("**Declarando Longitud y tiempo en el nuevo fichero**\n");
  if((retval=nc_def_dim(newNcId, "time", TIME, &newTimeDim)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "time", NC_DOUBLE, 1, &newTimeDim, &newTimeId)))
    ERR(retval);
  if((retval=nc_def_dim(newNcId, "longitude",LON, &newLonDim)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "longitude", NC_DOUBLE, 1, &newLonDim, &newLonId)))
    ERR(retval);
   if((retval=nc_def_dim(newNcId, "latitude", LAT, &newLatDim)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "latitude", NC_DOUBLE, 1, &newLatDim, &newLatId)))
    ERR(retval);
   //Declaramos la variable precipitacion ahora que tenemos latitud, longitud y tiempo
  dimIds[0]=newLatDim;
  dimIds[1]=newLonDim;
  dimIds[2]=newTimeDim;
  if((retval=nc_def_var(newNcId, "anomaliaClimatica", NC_SHORT, 3, dimIds, &newAnomaliaClimaticaId)))
    ERR(retval);
  if((retval = nc_enddef(newNcId)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newTimeId, time)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLonId, lon)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLatId, lat)))
    ERR(retval);
  printf("**Fichero creado con éxito**\n");
  free(lon);
  free(lat);
#pragma omp parallel private(retval,i,j,k,startp, countp, part_tg_val, estres_hidrico, estres_termico,anomalia_climatica, part_prec_val, device_drought_code, drought_code,last_drought_code,despl_lat,despl_total,timeId, orderedDroughtCode,orderedTg, precId,tgId, threadId, err, writeStartp, writeCountp, invalidValue) shared(ncIdPrec,ncIdTg,time,streams, device_part_tg_val, d_keys, newAnomaliaClimaticaId, newNcId,notCalculated)  num_threads(NUM_THREADS)
  {
    threadId = omp_get_thread_num();

    // Construct a new `nosync` execution policy with the custom stream
    auto nosync_exec_policy = thrust::cuda::par_nosync.on(streams[threadId]);
    if((retval=nc_inq_varid(ncIdPrec[threadId], "rr", &precId))) //Cogemos el id de la variable precipitacion
      ERR(retval);
    if((retval=nc_inq_varid(ncIdTg[threadId], "tg", &tgId))) //Cogemos el id de la variable temperatura
      ERR(retval);
    if((retval=nc_inq_varid(ncIdPrec[threadId], "time", &timeId))) //Cogemos el id de la variable tiempo del fichero de prec, aunque debería coincidir con el del fichero de tg.
      ERR(retval);
    printf("**Reservando memoria** \n");
    if((orderedDroughtCode = (float *)malloc(TIME*sizeof(float)*LON))==NULL)
      MEMERR("Reserva orderedPrec");
    if((orderedTg = (short *)malloc(TIME*sizeof(short)*LON)) == NULL)
      MEMERR("Reserva orderedTg");
    if((part_tg_val=(short *)calloc(sizeof(short), TIME*LON)) ==NULL)
      MEMERR("Reserva part_tg_val");
    //cudaMalloc((void**) &device_drought_code, TIME * sizeof(float));
    if((part_prec_val=(float *)calloc(sizeof(float),TIME*LON)) ==NULL)
      MEMERR("Reserva part_prec_val");
    if((drought_code = (float *)calloc(sizeof(float),TIME*LON))==NULL)
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
      *writeStartp = j; //Empezaremos a escribir en el fichero por la latitud que vayamos
      if(j == LAT/4)
	printf("**25%% del calculo completado**\n");
      else if(j == LAT/2)
	printf("**50%% del calculo completado**\n");
      else if(j == ((LAT/4)*3))
	printf("**75%% del calculo completado**\n");
      despl_lat = j * LON; //Desplazamiento Filas
      despl_total = despl_lat + k; //Desplazamiento total con desplazamiento de columna
      //Obtención prec y tg del lugar de todos los días.
      *startp = j;
      *(startp+1)= 0;//Siempre empezamos en la primera LON de la lat
      if((retval=nc_get_vara_short(ncIdTg[threadId],tgId,startp,countp, part_tg_val))){
	printf("error get tg\n");
	ERR(retval);
      }
      if((retval=nc_get_vara_float(ncIdPrec[threadId],precId,startp,countp, part_prec_val))){
	printf("error get prec\n");
	ERR(retval);
      }
      for(k=0;k<LON;k++){
	invalidValue[k] = calculateDroughtCode(part_tg_val+k*TIME, part_prec_val+k*TIME, drought_code+k*TIME, time);
      }
      cudaMemcpy(device_part_tg_val[threadId], part_tg_val,LON * TIME * sizeof(short),
		 cudaMemcpyHostToDevice);
      // --- Generate row indices: Las genera creando indices desde 0 hasta SizeToOrder * NUM_ITER, después crea un iterador con el SizeToOrder constante, y divide, de forma que queda solo el NUM_ITER correspondiente a cada elemento, es decir, la posición a la que pertenece e inicializa en d_keys
      thrust::transform(nosync_exec_policy,thrust::make_counting_iterator(0),
			thrust::make_counting_iterator(LON * TIME),
			thrust::make_constant_iterator(LON),
			d_keys[threadId],
			thrust::divides<short>());

      // --- Back-to-back approach
      //Ordenamos todos los elementos del vector
      thrust::stable_sort_by_key(nosync_exec_policy, device_part_tg_val[threadId],
				 device_part_tg_val[threadId] + LON * TIME,
				 d_keys[threadId],
				 thrust::less<short>());
      //Ordenamos las keys. De esta forma, al estar ordenado todo el vector, al ordenar las keys, cada posición estará ordenada
      // -6 2 4 1 3 5 ==> -6 1 2 3 4 5 ==> -6 2 4 1 3 5
      //  0 0 0 1 1 1 ==>  0 1 0 1 0 1 ==>  0 0 0 1 1 1
      thrust::stable_sort_by_key(nosync_exec_policy,d_keys[threadId],
				 d_keys[threadId]+LON*TIME,
				 device_part_tg_val[threadId],
				 thrust::less<short>());
    
    
      //Cálculo del estrés hídrico de cada día
      //Ordenamos los datos para calcular los percentiles
      memcpy(orderedDroughtCode, drought_code, TIME);
      for(i=0;i<LON;i++){
	if(!invalidValue[i]){
	  thrust::sort(thrust::host, orderedDroughtCode+i*TIME, orderedDroughtCode+i*TIME + TIME);
	  calculateEstresHidrico(drought_code+i*TIME,orderedDroughtCode+i*TIME, estres_hidrico+i*TIME);
	}
      }
      //Ordenamos para el cálculo de los percentiles
      cudaMemcpy(orderedTg, device_part_tg_val[threadId], TIME * sizeof(short) * LON,
		 cudaMemcpyDeviceToHost);
      //Cálculo del estrés térmico de cada día
      for(i=0;i<LON;i++){
	if(!invalidValue[i]){
	  calculateEstresTermico(part_tg_val+i*TIME, orderedTg+i*TIME, estres_termico+i*TIME);
	  //Cálculo de la anomalía climática
	  calculateAnomaliaClimatica(estres_termico+i*TIME, estres_hidrico+i*TIME, anomalia_climatica+i*TIME);
	} else{
	  printf("Invalid values in LAT: %d LON: %d",j,k);
	  #pragma omp atomic
	  notCalculated++;
	  for(k=0; k<TIME; k++){
	    *(estres_termico+i*TIME+k)=-1;
	    *(estres_hidrico+i*TIME+k)=-1;
	    *(anomalia_climatica+i*TIME+k) = -1;
	  }
	}
      }
      //Escribimos en el nuevo fichero la anomalia climatica
      *(writeCountp+2)=TIME; // Escribimos todos los datos de un dia
      *(writeCountp+1)=LON; //Escribimos un dato de LON
      *(writeCountp)=1; //Escribimos un dato de LAT
      *(writeStartp+1)=0;//Empezamos a escribir desde la LON que toque
      *(writeStartp+2)=0;//Siempre empezamos a escribir desde el primer dia
#pragma omp critical
      {
	nc_put_vara_short(newNcId, newAnomaliaClimaticaId,writeStartp, writeCountp, anomalia_climatica);
      }
    }
    printf("**Cálculos finalizados con éxito proceso: %d** \n", omp_get_thread_num());
    cudaStreamDestroy(streams[threadId]);
    free(orderedTg);
    free(orderedDroughtCode);
    free(startp);
    free(countp);
    free(drought_code);
    free(part_prec_val);
    free(part_tg_val);
    free(estres_termico);
    free(estres_hidrico);
    free(anomalia_climatica);
    cudaFree(device_part_tg_val);
    cudaFree(d_keys);
    //   cudaFree(device_drought_code);
    nc_close(ncIdPrec[threadId]);
    nc_close(ncIdTg[threadId]);
  }
  
  gettimeofday(&t2, NULL);
  segundos = (((t2.tv_usec - t.tv_usec)/1000000.0f)  + (t2.tv_sec - t.tv_sec));
 
  printf("\n *******> Duracion total: %f segundos\n\n", segundos);
   printf("Esperando a cerrar el archivo creado\n");
  if((retval = nc_close(newNcId)))
    ERR(retval);
  printf("Total:%d , NotCalculated: %ld\n",LAT*LON, notCalculated);
  return 0;
}

int calculateDroughtCode(short *temperature, float *rainfall, float *droughtCode, double *days){
  int i, limite=(LON*TIME) - ((LON*TIME/100)*70); //Si un 30% de los datos se ha tenido que suponer, ignorar esta posición
  int invalid, totalInvalid=0;
  float prevDroughtCode= 0;
  float efRainfall, moisture, prevMoisture, evapotranspiration;
  for(i=1;i<TIME;i++){
    invalid=0;
    if(*(rainfall+i) < 0.0 || *(rainfall+i) > 300.0){
      *(rainfall+i)=0; //Si no se obtienen datos, se asume que no ha llovido
      invalid=1;
    }
    if(*(temperature+i) < -8000 || *(temperature+i) > 8000){
      *(temperature+i) = 1040; //Temperatura media del dataset
      invalid=1;
    }
    if(invalid)
      totalInvalid++;
    if(*(rainfall+i)>2.8){
      efRainfall = 0.86 * *(rainfall+i) -1.27;
      prevMoisture = 800 * exp(-*(droughtCode+i)/400);
      moisture = prevMoisture + 3.937 * efRainfall;
      prevDroughtCode = 400 * log(800/moisture);
      prevDroughtCode = prevDroughtCode < 0 ? 0 : prevDroughtCode;
    }
    evapotranspiration = 0.36 * (*temperature + 2.0) + getDayLength(*(days+i));
    evapotranspiration = evapotranspiration < 0 ? 0 : evapotranspiration;
    *(droughtCode+i) = prevDroughtCode + 0.5 * evapotranspiration;
  }
  if(limite < totalInvalid)
    return 1;
  return 0;
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
    *(estres_termico+i) =-1;
  for(i=7; i<TIME; i++){
    media_ocho_dias = (short) (((int) *(part_tg_val+i-7) + (int) *(part_tg_val+i-6) + (int) *(part_tg_val+i-5) + (int) *(part_tg_val+i-4) + (int) *(part_tg_val+i-3) + (int) *(part_tg_val+i-2) + (int) *(part_tg_val+i-1) + (int) *(part_tg_val+i)) /8);
    *(estres_termico+i) = encontrarPercentilShort(percentiles,media_ocho_dias);
  }
  free(percentiles);
  return 0;
}

int calculateEstresHidrico(float *part_prec_val,float *orderedPrec, short *estres_hidrico){
  int i;
  int tamano_percentil = TIME/100;
  float *percentiles = (float *)malloc(99*sizeof(float));
  //Sacamos los números que dividen los percentiles
  for(i=1; i<99; i++)
    percentiles[i-1]=*(orderedPrec+i*tamano_percentil);
  //Sacamos los números que dividen los percentiles
  for(i=1; i<99; i++)
    percentiles[i-1]=*(orderedPrec+i*tamano_percentil);
  for(i=0; i<TIME; i++){ //Se necesita la media de los 8 días anteriores
    *(estres_hidrico+i) = encontrarPercentilFloat(percentiles, *(part_prec_val+i));
  }
  free(percentiles);
  return 0;
}

int calculateAnomaliaClimatica(short *estres_termico, short *estres_hidrico, short* anomalia_climatica){ 
  int i = 0;
  
  for(i=0; i<TIME; i++)
    if(*(estres_termico+i)<1 ||*(estres_termico+i)>99 || *(estres_hidrico+i)<1 ||*(estres_hidrico+i)>99)
      *(anomalia_climatica+i) = -1;
    else
      *(anomalia_climatica+i) = (*(estres_termico+i) + *(estres_hidrico+i)) / 2;
  return 0;
}

short encontrarPercentilShort(short *percentiles, short val){
  short last_percentil=0, first=0, middle=49, last=99;
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
short encontrarPercentilFloat(float *percentiles, float val){
  short last_percentil=0, first=0, middle=49, last=99;
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
