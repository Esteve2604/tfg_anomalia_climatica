#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
//Definitions from netCDF library
#define NC_SHORT        3
#define NC_DOUBLE       6
#define NC_FLOAT        5

#define FILE_NAME_ANOMALIA "./copernicus_anomalia_climatica_1995_2010.nc"
#define ERRCODE 2
#define MEMERRORCODE 3
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define MEMERR(e) {printf("Error: %s\n", e); exit(MEMERRORCODE);}

#define NUM_ESCR_ANOMALIA 1
#define NDIMS 3
#define LAT 437
#define LON 592
#define TIME 5844

int write_anomalia(int ncIdAnomalia, int anomaliaId, int newNcId, int newAnomaliaId);

int main(){
  int i,ncIdAnomalia, ncIdTg, lonId, latId, timeId, anomaliaId, retval;
  int newNcId, newLonDim, newLatDim, newTimeDim, dimIds[3], newAnomaliaId, newLonId, newLatId, newTimeId; //Ida for new file
  double lat[LAT];
  double lon[LON];
  double time[TIME];
  struct timeval t, t2;
  double segundos;
  gettimeofday(&t, NULL);

  //Primero vamos a optimizar la precicpitación
  //Lectura de datos
  printf("**Creando el fichero limpio de anomalia**\n");
  if((retval=nc_create("copernicus_anomalia_climatica_1995_2010_by_days.nc", NC_CLOBBER, &newNcId)))
    ERR(retval);
  printf("**Fichero creado con éxito**\n");
  printf("**Lectura de datos**\n");
  if((retval=nc_open(FILE_NAME_ANOMALIA, NC_NOWRITE, &ncIdAnomalia)))
    ERR(retval);
  if((retval=nc_inq_varid(ncIdAnomalia, "longitude", &lonId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdAnomalia, lonId, lon)))
    ERR(retval);
  if((retval=  nc_inq_varid(ncIdAnomalia, "latitude", &latId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdAnomalia, latId, lat)))
    ERR(retval);
  if((retval= nc_inq_varid(ncIdAnomalia, "time", &timeId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdAnomalia, timeId, time)))
    ERR(retval);
  if((retval=nc_inq_varid(ncIdAnomalia, "anomaliaClimatica", &anomaliaId)))
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
  printf("**Declarando Latitud  y prec en el nuevo fichero, y escribiendo Lat,lon y time**\n");
  if((retval=nc_def_dim(newNcId, "latitude", LAT, &newLatDim)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "latitude", NC_DOUBLE, 1, &newLatDim, &newLatId)))
    ERR(retval);
  //Declaramos la variable precipitacion ahora que tenemos latitud, longitud y tiempo
  dimIds[1]=newLatDim;
  dimIds[2]=newLonDim;
  dimIds[0]=newTimeDim;
  if((retval=nc_def_var(newNcId, "anomaliaClimatica", NC_SHORT, 3, dimIds, &newAnomaliaId)))
    ERR(retval);
  if((retval = nc_enddef(newNcId)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newTimeId, time)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLonId, lon)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLatId, lat)))
    ERR(retval);
  //NO Liberamos la latitud y longitud nuevas para ponerlas en el fichero de temperatura
  for(i=0;i<LAT;i++)
    printf("%f\t",lat[i]);
  printf("**Escritura Correcta**\n");
  write_anomalia(ncIdAnomalia, anomaliaId, newNcId, newAnomaliaId);
  if((retval = nc_close(newNcId)))
    ERR(retval);
  printf("**Escritura de anomalia en fichero nuevo completado**\n");
  printf("**SUCCESS!!**\n");
  gettimeofday(&t2, NULL);
  segundos = (((t2.tv_usec - t.tv_usec)/1000000.0f)  + (t2.tv_sec - t.tv_sec));
  printf("\n *******> Duracion total: %f segundos\n\n", segundos);
  return 0;
}

/* Debido a que se puede llegar a leer un gran volumen de datos, y vamos a transformar los datos de forma que esten de forma optimizada
para acceder de forma temporal, se ha decidido dividir las escrituras. Como las escrituras deben ser directas, se ha establecido una variable
global que establece en cuantas veces se dividirán las escrituras. Divide el número de latitudes que se escribirán a la vez. NUM_ESCR_PREC
*/
int write_anomalia(int ncIdAnomalia, int anomaliaId, int newNcId, int newAnomaliaId){
 int i,j,k,z,writeDayIndex=0;
  int retval;
  unsigned long newAnomaliaIndex =0, despl_total;
  size_t *startp,*countp, *writeStartp, *writeCountp;
  short *partAnomaliaVal;
  short *escrAnomaliaVal;
  int lastDayThisIter=0;
  int numDaysPerIter= TIME/NUM_ESCR_ANOMALIA;
  //Vamos a leer 1 Lat a la vez para poder interpolar con los datos de alrededor.
  if((partAnomaliaVal=malloc(sizeof(short)*TIME)) == NULL)
    MEMERR("Reserva de tgVal para LAT");
  //Reservamos necesario para el máximo número posible ha escribir
  if((escrAnomaliaVal=malloc(sizeof(short)*((numDaysPerIter+TIME%NUM_ESCR_ANOMALIA)*LAT*LON))) == NULL)
    MEMERR("Reserva de escrTgVal para LAT");
  
  //Para lecturas. Orden del fichero de lectura  LAT -> LON -> TIME
  startp=(size_t *)malloc(NDIMS*sizeof(size_t));
  countp=(size_t *)malloc(NDIMS*sizeof(size_t));
  *countp=1; //Leemos un dato de LAT
  *(countp+1)=1; // Un dato de LON
  *(countp+2)=TIME; //todos los datos de TIME
  
  //Para escrituras. Orden del fichero nuevo TIME -> LAT -> LON
  writeStartp=(size_t *)malloc(NDIMS*sizeof(size_t));
  writeCountp=(size_t *)malloc(NDIMS*sizeof(size_t));
  *(writeCountp+2)=LON; // Escribimos todos los datos de un LON
  *(writeCountp+1)=LAT; //Escribimos todos los datos de LAT
  *(writeStartp+1)=0;//Siempre empezamos a escribir desde la primera LAT
  *(writeStartp+2)=0;//Siempre empezamos a escribir desde el primer LON
  printf("El fichero se ha dividido en %d para escribir\n",NUM_ESCR_ANOMALIA);
  for(z=0; z<NUM_ESCR_ANOMALIA; z++){
    lastDayThisIter=numDaysPerIter*(z+1);//La ultima LAT de esta iter
    if(z==NUM_ESCR_ANOMALIA-1){//Si es la ultima hacemos las sobrantes
      lastDayThisIter=TIME;
    }
    for(i=0; i<LAT; i++){
      *startp=i;
      //Marcamos las latitudes que nos tocan esta iter
      for(j=0; j<LON; j++){
	*(startp+1)=j;
	//Leemos la anomalia.
	if((retval=nc_get_vara_short(ncIdAnomalia,anomaliaId,startp,countp, partAnomaliaVal))){
	  printf("error get anomalia \n");
	  ERR(retval);
	}
	newAnomaliaIndex=i*LON+j;//Iniciamos el indice por la LAT por la que vamos
	for(k=z*numDaysPerIter; k<lastDayThisIter; k++){
	  //Desplazamos el puntero para que todas las posiciones tengan sus días contiguos.
	  *(escrAnomaliaVal+newAnomaliaIndex)=  *(partAnomaliaVal+k);
	  newAnomaliaIndex = newAnomaliaIndex + LAT*LON;
	}	
      }
    }
    //El indice por el que tenemos que empezar a insertar
    *writeStartp=writeDayIndex;
    //El numero de dias que tenemos que escribir
    *writeCountp=lastDayThisIter - z*numDaysPerIter;
    //Escribimos en el fichero
    nc_put_vara_short(newNcId, newAnomaliaId,writeStartp, writeCountp, escrAnomaliaVal);
    //Actualizamos el indice por el que empezaremos a escribir la proxima iter
    writeDayIndex=lastDayThisIter;
    printf("%d de %d partes del fichero escritas\n",z+1, NUM_ESCR_ANOMALIA);
  }
  free(partAnomaliaVal);
  free(escrAnomaliaVal);
  free(startp);
  free(countp);
  free(writeStartp);
  free(writeCountp);
  return 0;
}
