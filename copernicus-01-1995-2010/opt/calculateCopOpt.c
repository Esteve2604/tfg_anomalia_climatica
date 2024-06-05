#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#define FILE_NAME_PREC "../copernicus_prec_data_clean.nc"
#define FILE_NAME_TG "../copernicus_tg_data_clean.nc"
#define NDIMS 3
#define LAT 437
#define LON 592
#define TIME 5844

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
int calculateEstresTermico(short *part_tg_val, short *estres_termico);

//Deja los datos calculados en estres_hidrico y ordena de menor a mayor part_prec_val
int calculateEstresHidrico(double *part_prec_val, short *estres_hidrico);

//Deja los datos calculados en anomalia_climatica
int calculateAnomaliaClimatica(short *estres_termico, short *estres_hidrico, short* anomalia_climatica);
//Función para ordenar shorts con qsort
static int shortCompare(const void *p1, const void *p2);
//Función para ordenar double con qsort
static int doubleCompare(const void *p1, const void *p2);

//Función para encontrar el percentil al que pertenece un valor, devuelve el percentil
short encontrarPercentilShort(short *percentiles, short val);
short encontrarPercentilDouble(double *percentiles, double val);

int main(){
  int retval;
  int i,j,k, ncIdPrec, ncIdTg;
  size_t *startp, *countp;
  float varValue;
  short *part_tg_val, *estres_hidrico, *estres_termico, *anomalia_climatica;
  double *part_prec_val, *drought_code, last_drought_code;
  int despl_time, despl_lat, despl_total;
  int precId, tgId, timeId;
  int invalidValue, notCalculated;
  double time[TIME];
  struct timeval t, t2;
  double segundos;
  gettimeofday(&t, NULL);
  printf("**Inicio** \n");
  printf("**Abriendo ficheros** \n");
  if((retval= nc_open(FILE_NAME_PREC, NC_NOWRITE, &ncIdPrec)))
    ERR(retval);
  if((retval=nc_open(FILE_NAME_TG, NC_NOWRITE, &ncIdTg)))
    ERR(retval);
  if((retval=nc_inq_varid(ncIdPrec, "rr", &precId)))
    ERR(retval);
  if((retval=nc_inq_varid(ncIdTg, "tg", &tgId)))
    ERR(retval);
  if((retval=nc_inq_varid(ncIdPrec, "time", &timeId)))
    ERR(retval);
  printf("**Reservando memoria** \n");
  if((part_tg_val=calloc(sizeof(short), TIME)) ==NULL)
    MEMERR("Reserva part_tg_val");
  if((part_prec_val=calloc(sizeof(double),TIME)) ==NULL)
    MEMERR("Reserva part_prec_val");
  if((drought_code = calloc(sizeof(double),TIME))==NULL)
    MEMERR("Reserva drought_code");
  if((estres_termico = calloc(sizeof(short),TIME))==NULL)
    MEMERR("Reserva estres_termico");
  if((estres_hidrico = calloc(sizeof(short),TIME))==NULL)
    MEMERR("Reserva estres_hidrico");
  if((anomalia_climatica = calloc(sizeof(short),TIME))==NULL)
    MEMERR("Reserva anomalia_climatica");
  startp=malloc(NDIMS*sizeof(size_t));
  countp=malloc(NDIMS*sizeof(size_t));
  *countp=1; //Leemos un dato de lat
  *(countp+1)=1; //un dato de lon
  *(countp+2)=TIME; //y tantos datos como dias
  *(startp+2)=0;//Siempre empezamos a leer desde el primer dia
  nc_get_var_double(ncIdPrec, timeId, time);
  i=0;
  for(j=0;j<LAT; j++){
    if(j == LAT/4)
      printf("**25%% del calculo completado**\n");
    else if(j == LAT/2)
      printf("**50%% del calculo completado**\n");
    else if(j == ((LAT/4)*3))
      printf("**75%% del calculo completado**\n");
    despl_lat = j * LON; //Desplazamiento Filas
    for(k=0;k<LON;k++){
      invalidValue=0;
      despl_total = despl_lat + k; //Desplazamiento total con desplazamiento de columna
      //Obtención prec y tg del lugar de todos los días.
      *startp = j;
      *(startp+1)= k;
      if((retval=nc_get_vara_double(ncIdPrec,precId,startp,countp, part_prec_val)))
	ERR(retval);
      if((retval=nc_get_vara_short(ncIdTg,tgId,startp,countp, part_tg_val)))
	ERR(retval);
      //Calculo del drought code
      last_drought_code=calculateDroughtCode(part_tg_val[0], part_prec_val[0], 0, time[0]);
      if(last_drought_code<0)
	invalidValue = 1;
      else {
	drought_code[0]=last_drought_code;
	for(i=1;i<TIME;i++){
	  last_drought_code=calculateDroughtCode(part_tg_val[i], part_prec_val[i], last_drought_code,time[i]);
	  drought_code[i]=last_drought_code;
	  if(last_drought_code<0){
	    invalidValue = 1;
	    break;
	  }		
	}
      }
      if(!invalidValue){
	//Cálculo del estrés térmico de cada día
	calculateEstresTermico(part_tg_val, estres_termico);
     
	//Cálculo del estrés hídrico de cada día
	calculateEstresHidrico(drought_code, estres_hidrico);
      
	//Cálculo de la anomalía climática
	calculateAnomaliaClimatica(estres_termico, estres_hidrico, anomalia_climatica);
      } else {
	notCalculated++;
        for(i=0; i<TIME; i++){
	  *(estres_termico+i)=-1;
	  *(estres_hidrico+i)=-1;
	  *(anomalia_climatica+i) = -1;
	}
      }
    }
  }
  printf("**Cálculos finalizados con éxito** \n");
  gettimeofday(&t2, NULL);
  segundos = (((t2.tv_usec - t.tv_usec)/1000000.0f)  + (t2.tv_sec - t.tv_sec));

  printf("\n *******> Duracion total: %f segundos\n\n", segundos);
  free(startp);
  free(countp);
  free(drought_code);
  free(part_prec_val);
  free(part_tg_val);
  free(estres_termico);
  free(estres_hidrico);
  free(anomalia_climatica);
  nc_close(ncIdPrec);
  nc_close(ncIdTg);
  printf("Total:%d , NotCalculated: %d\n",LAT*LON, notCalculated);
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
}

int calculateEstresTermico(short *part_tg_val, short *estres_termico){
  int i;
  int tamano_percentil = TIME/100;
  short *percentiles = malloc(99*sizeof(short));
  short *orderedTg = malloc(TIME*sizeof(part_tg_val));
  short media_ocho_dias;
  memcpy(orderedTg, part_tg_val, TIME);
  //Ordenamos los datos para calcular los percentiles
  qsort(orderedTg,TIME, sizeof(short), shortCompare);
  //Sacamos los números que dividen los percentiles
  for(i=1; i<99; i++)
    percentiles[i-1]=*(orderedTg+i*tamano_percentil);
  for(i=0; i<7; i++) //Se necesita la media de los 8 días anteriores
    *(estres_termico) =-1;
  for(i=7; i<TIME; i++){
    media_ocho_dias = (short) (((int) *(part_tg_val+i-7) + (int) *(part_tg_val+i-6) + (int) *(part_tg_val+i-5) + (int) *(part_tg_val+i-4) + (int) *(part_tg_val+i-3) + (int) *(part_tg_val+i-2) + (int) *(part_tg_val+i-1) + (int) *(part_tg_val+i)) /8);
    *(estres_termico) = encontrarPercentilShort(percentiles, media_ocho_dias);
  }
  free(percentiles);
  free(orderedTg);
  return 0;
}

int calculateEstresHidrico(double *part_prec_val, short *estres_hidrico){
  int i;
  int tamano_percentil = TIME/100;
  double *percentiles = malloc(99*sizeof(double));
  double *orderedPrec = malloc(TIME*sizeof(part_prec_val));
  memcpy(orderedPrec, part_prec_val, TIME);
  //Ordenamos los datos para calcular los percentiles
  qsort(orderedPrec,TIME, sizeof(double), doubleCompare);
  //Sacamos los números que dividen los percentiles
  for(i=1; i<99; i++)
    percentiles[i-1]=*(orderedPrec+i*tamano_percentil);
  for(i=0; i<TIME; i++){
    *(estres_hidrico) = encontrarPercentilDouble(percentiles, *(part_prec_val+i));
  }
  free(percentiles);
  free(orderedPrec);
  return 0;
}

int calculateAnomaliaClimatica(short *estres_termico, short *estres_hidrico, short* anomalia_climatica){
  int i = 0;
  for(i=0; i<TIME; i++)
    *(anomalia_climatica+i) = (*(estres_termico+i) + *(estres_hidrico+i)) / 2;
  return 0;
}

//Función para ordenar shorts con qsort
static int shortCompare(const void *p1, const void *p2){
  short short_a = *((short *)p1);
  short short_b = *((short *)p2);
  if(short_a == short_b)
    return 0;
  else if(short_a < short_b)
    return -1;
  else
    return 1;
}
//Función para ordenar double con qsort
static int doubleCompare(const void *p1, const void *p2){
  double double_a = *((double *)p1);
  double double_b = *((double *)p2);
  if(double_a == double_b)
    return 0;
  else if(double_a < double_b)
    return -1;
  else
    return 1;
}

short encontrarPercentilShort(short *percentiles, short val){
  short last_percentil=0, first=0, middle=49, last=98;
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
short encontrarPercentilDouble(double *percentiles, double val){
  short last_percentil=0, first=0, middle=49, last=98;
  while(first<=last){
    if(*(percentiles+middle) < val){
      last_percentil=middle; //Guardamos el percentil inferior
      first = middle +1;
    } else if(*(percentiles+middle)== val){
      last_percentil = middle; 
      break;
    } else {
      last = middle -1;
    }
    middle = (first +last) /2;
  }
  return last_percentil +1; //+1 porque la posición del array empieza en 0
}
