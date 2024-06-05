#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

//Definitions from netCDF library
#define NC_SHORT        3
#define NC_DOUBLE       6


#define FILE_NAME_PREC "./rr_ens_mean_0.1deg_reg_v29.0e.nc"
#define FILE_NAME_TG "./tg_ens_mean_0.1deg_reg_v29.0e.nc"

#define ERRCODE 2
#define MEMERRORCODE 3
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define MEMERR(e) {printf("Error: %s\n", e); exit(MEMERRORCODE);}

#define NDIMS 3
#define LAT 465
#define LON 705
#define TIME 27028

int check_lon(double *lon, int lonSize, int latSize, int timeSize, double *val, double *finalLon, int *validLon); //Retorna en valid el número de valores validos final
int check_lat(double *lat, int lonSize, int latSize , int timeSize, double *val,double *finalLat, int *validLat);
int check_prec(double *prec, double *oldLon, double *oldLat, int lonSize, int latSize, int timeSize, double *newPrec);
int check_tg(short *tg, double *oldLon, double *oldLat, int lonSize, int latSize, int timeSize, short *newTg);
int main(){
  int retval;
  int ncIdPrec, ncIdTg , lonId, latId,timeId, precId, tgId, newLonSize, newLatSize;
  int newNcId, newLonDim, newLatDim, newTimeDim, dimIds[3], newPrecId, newTgId, newLonId, newLatId, newTimeId; //Id for new file
  double *precVal, *newLon, *newLat, *newPrec ;
  short *tgVal, *newTg;
  double lat[LAT];
  double lon[LON];
  double time[TIME];
  struct timeval t, t2;
  double segundos;
  gettimeofday(&t, NULL);
  //Lectura de datos
  printf("**Lectura de datos**\n");
  if((retval=nc_open(FILE_NAME_PREC, NC_NOWRITE, &ncIdPrec)))
    ERR(retval);
  if((retval=nc_open(FILE_NAME_TG, NC_NOWRITE, &ncIdTg)))
    ERR(retval);
  if((retval=nc_inq_varid(ncIdPrec, "longitude", &lonId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdPrec, lonId, lon)))
    ERR(retval);
  if((retval=  nc_inq_varid(ncIdPrec, "latitude", &latId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdPrec, latId, lat)))
    ERR(retval);
  if((retval= nc_inq_varid(ncIdPrec, "time", &timeId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdPrec, timeId, time)))
    ERR(retval);
  if((precVal = calloc(sizeof(double), LON * LAT * TIME)) == NULL)
    MEMERR("Reserva de precVal, linea 54");
  if((tgVal = calloc(sizeof(short), LON * LAT * TIME)) == NULL)
    MEMERR("Reserva de tgVal, linea 56");
 
  if((retval=nc_inq_varid(ncIdPrec, "rr", &precId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdPrec, precId, precVal)))
    ERR(retval);
  if((retval=nc_inq_varid(ncIdTg, "tg", &tgId)))
    ERR(retval);
  if((retval=nc_get_var_short(ncIdTg, tgId, tgVal)))
    ERR(retval);
  printf("**Lectura de datos Correcta**\n");
  //Comprobación y corrección de datos
  printf("**Comprobando la Longitud**\n");
  if((newLon = malloc(LON * sizeof(double))) == NULL)
    MEMERR("Reserva de newLon, linea 70");
  check_lon(lon, LON,LAT,TIME, precVal, newLon, &newLonSize);
  printf("**Comprobando la Latitud**\n");
  if((newLat = malloc(LAT * sizeof(double))) == NULL)
    MEMERR("Reserva de newLat, linea 75");
  check_lat(lat, LON,LAT,TIME, precVal, newLat, &newLatSize);  
  printf("**Comprobando la Precipitacion**\n");
  if((newPrec = malloc(newLonSize * newLatSize * TIME *sizeof(double))) == NULL)
    MEMERR("Reserva de newPrec, linea 71");
  check_prec(precVal,lon, lat,LON,LAT,TIME, newPrec);
  free(precVal);//Liberamos la memoria de los datos ya procesados
  printf("**Comprobando la temperatura**\n");
  if((newTg = malloc(newLonSize * newLatSize * TIME *sizeof(short))) == NULL)
    MEMERR("Reserva de newTg, linea 75");
  check_tg(tgVal,lon, lat, LON, LAT, TIME, newTg);
  free(tgVal); //Liberamos la memoria de los datos ya procesados
  printf("**Todas las comprobaciones han sido correctas**\n");
  //Creando el nuevo fichero
  //Precipitacion
  printf("**Creando el fichero limpio de precipitacion**\n");
  if((retval=nc_create("copernicus_prec_data_clean.nc", NC_CLOBBER, &newNcId)))
    ERR(retval);
  if((retval=nc_def_dim(newNcId, "time", TIME, &newTimeDim)))
    ERR(retval);
  if((retval=nc_def_dim(newNcId, "longitude",newLonSize, &newLonDim)))
    ERR(retval);
  if((retval=nc_def_dim(newNcId, "latitude", newLatSize, &newLatDim)))
    ERR(retval);
  
  dimIds[0]=newLatDim;
  dimIds[1]=newLonDim;
  dimIds[2]=newTimeDim;
  if((retval=nc_def_var(newNcId, "time", NC_DOUBLE, 1, &newTimeDim, &newTimeId)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "longitude", NC_DOUBLE, 1, &newLonDim, &newLonId)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "latitude", NC_DOUBLE, 1, &newLatDim, &newLatId)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "rr", NC_DOUBLE, 3, dimIds, &newPrecId)))
    ERR(retval);
  if((retval = nc_enddef(newNcId)))
    ERR(retval);
  //Escritura de datos
  printf("**Escribiendo datos**\n");
  if ((retval = nc_put_var_double(newNcId, newTimeId, time)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLonId, newLon)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLatId, newLat)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newPrecId, newPrec)))
    ERR(retval);
 
  //Cierre
  if((retval = nc_close(newNcId)))
    ERR(retval);
  printf("**Creacion y llenado de precipitacion completado**\n");
  //Temperatura
  //Creación de fichero
  printf("**Creando el fichero limpio de temperatura**\n");
  if((retval=nc_create("copernicus_tg_data_clean.nc", NC_CLOBBER, &newNcId)))
    ERR(retval);
  if((retval=nc_def_dim(newNcId, "time", TIME, &newTimeDim)))
    ERR(retval);
  if((retval=nc_def_dim(newNcId, "longitude",newLonSize, &newLonDim)))
    ERR(retval);
  if((retval=nc_def_dim(newNcId, "latitude", newLatSize, &newLatDim)))
    ERR(retval);
  
  dimIds[0]=newLatDim;
  dimIds[1]=newLonDim;
  dimIds[2]=newTimeDim;
  if((retval=nc_def_var(newNcId, "longitude", NC_DOUBLE, 1, &newLonDim, &newLonId)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "latitude", NC_DOUBLE, 1, &newLatDim, &newLatId)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "time", NC_DOUBLE, 1, &newTimeDim, &newTimeId)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "tg", NC_SHORT, 3, dimIds, &newTgId)))
    ERR(retval);
  if((retval = nc_enddef(newNcId)))
    ERR(retval);
  //Escritura de Datos
  printf("**Escribiendo datos**\n");
  if ((retval = nc_put_var_double(newNcId, newTimeId, time)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLonId, newLon)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLatId, newLat)))
    ERR(retval);
  if ((retval = nc_put_var_short(newNcId, newTgId, newTg)))
    ERR(retval);
 
  //Cierre
  if((retval = nc_close(newNcId)))
    ERR(retval);
  
  free(newLon);
  free(newLat);
  free(newPrec);
  free(newTg);
  
  printf("**SUCCESS!!**\n");
  gettimeofday(&t2, NULL);
  segundos = (((t2.tv_usec - t.tv_usec)/1000000.0f)  + (t2.tv_sec - t.tv_sec));
  
  printf("\n *******> Duracion total: %f segundos\n\n", segundos);
  return 0;
}

int check_lon(double *lon,int lonSize, int latSize, int timeSize, double *val, double *finalLon, int *lonValid){
  int i,j,k;
  unsigned long despl_part,despl_total;
  unsigned long limite =(latSize*timeSize) - ((latSize*timeSize/100)*5);
  unsigned long notValid;
  *lonValid=0;
  for(i=0; i<lonSize; i++){
    notValid=0;
    for(j=0; j<latSize; j++){
      despl_part= j*lonSize + i; //Desplazamiento de fila y columna
      for(k=0; k<timeSize; k++){
	despl_total = k*latSize*lonSize+despl_part; //"Cambio matriz de dia"
	if(*(val+despl_total) < 0.0 || *(val+despl_total) > 300.0)
	  notValid++;
      }
    }
    if(notValid < limite){
      *(finalLon+(*lonValid)) = *(lon+i);
      *(lonValid) = *(lonValid) +1;
    } else {
      *(lon+i) = -9999.0;
    }
  }
  return 0;
}
int check_lat(double *lat,int lonSize, int latSize, int timeSize, double *val, double *finalLat, int *latValid){
  int i,j,k;
  unsigned long despl_part,despl_total;
  unsigned long limite =(lonSize*timeSize) - ((lonSize*timeSize/100)*5);
  unsigned long notValid;
  *latValid=0;
  for(i=0; i<latSize; i++){
    notValid=0;
    for(j=0; j<lonSize; j++){
      despl_part= i*lonSize + j; //Desplazamiento de fila y columna
      for(k=0; k<timeSize; k++){
	despl_total = k*latSize*lonSize+despl_part; //"Cambio de matriz de dia"
	if(*(val+despl_total) < 0 || *(val+despl_total) > 300.0)
	  notValid++;
      }
    }
    if(notValid < limite){
      *(finalLat+(*latValid)) = *(lat+i);
      *(latValid) = *(latValid) +1;
    } else{
      *(lat+i) = -9999.0;
    }
  }
  return 0;
}
int check_prec(double *prec,double *oldLon, double *oldLat, int lonSize, int latSize, int timeSize, double *newPrec){
  int i,j,k;
  unsigned long despl_part,despl_total, newPrecIndex =0;
  for(i=0; i<latSize; i++){
    if(*(oldLat+i) == -9999.0)
      continue;
    for(j=0; j<lonSize; j++){
      if(*(oldLon+j) == -9999.0)
	continue;
      despl_part= i*lonSize + j; //Desplazamiento de fila y columna
      for(k=0; k<timeSize; k++){
	despl_total = k*latSize*lonSize+despl_part; //"Cambio de matriz de dia"
	if(*(prec+despl_total) < 0){ //En caso de que sea nulo, copiamos el valor de un lugar cercano
	  if(*(prec+despl_total-1) >= 0.0) //A la izquierda
	    *(prec+despl_total) = *(prec+despl_total-1);
	  else if(*(prec+despl_total+1) >= 0.0) //Derecha
	    *(prec+despl_total) = *(prec+despl_total+1);
	  else if(*(prec+despl_total+lonSize) >= 0.0) //Abajo
	    *(prec+despl_total) = *(prec+despl_total+lonSize);
	  else if(*(prec+despl_total-lonSize) >= 0.0) //Arriba
	    *(prec+despl_total) = *(prec+despl_total+lonSize);
	}
	*(newPrec+newPrecIndex)= *(prec+despl_total);
	newPrecIndex++;
      }
    }
  }
  return 0;
}
int check_tg(short *tg, double *oldLon, double *oldLat, int lonSize, int latSize, int timeSize, short *newTg){
  int i,j,k, newTgIndex =0;
  unsigned long despl_part,despl_total, newPrecIndex =0;
  for(i=0; i<latSize; i++){
    if(*(oldLat+i) == -9999.0)
      continue;
    for(j=0; j<lonSize; j++){
      if(*(oldLon+j) == -9999.0)
	continue;
      despl_part= i*lonSize + j; //Desplazamiento de fila y columna
      for(k=0; k<timeSize; k++){
	despl_total = k*latSize*lonSize+despl_part; //"Cambio de matriz de dia"
	if(*(tg+despl_total) < -8500 || *(tg+despl_total) > 8500){ //En caso de que sea nulo, copiamos el valor de un lugar cercano
	  if(*(tg+despl_total-1) > -8500 && *(tg+despl_total-1) < 8500) //Arriba
	    *(tg+despl_total) = *(tg+despl_total-1);
	  else if(*(tg+despl_total+1) > -8500 && *(tg+despl_total+1) < 8500) //Abajo
	    *(tg+despl_total) = *(tg+despl_total+1);
	  else if(*(tg+despl_total+lonSize) > -8500 && *(tg+despl_total+lonSize) < 8500) //Derecha
	    *(tg+despl_total) = *(tg+despl_total+lonSize);
	  else if(*(tg+despl_total-lonSize) > -8500 && *(tg+despl_total-lonSize) < 8500) //Izquierda
	    *(tg+despl_total) = *(tg+despl_total-lonSize);
	}
	*(newTg+newTgIndex)= *(tg+despl_total);
	newTgIndex++;
      }
    }
  }
  return 0;
}
