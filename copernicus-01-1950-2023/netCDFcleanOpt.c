#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
//Definitions from netCDF library
#define NC_SHORT        3
#define NC_DOUBLE       6
#define NC_FLOAT        5

#define FILE_NAME_PREC "./rr_ens_mean_0.1deg_reg_v29.0e.nc"
#define FILE_NAME_TG "./tg_ens_mean_0.1deg_reg_v29.0e.nc"
#define ERRCODE 2
#define MEMERRORCODE 3
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define MEMERR(e) {printf("Error: %s\n", e); exit(MEMERRORCODE);}

#define NUM_ESCR_PREC 8
#define NUM_ESCR_TG 4
#define NDIMS 3
#define LAT 465
#define LON 705
#define TIME 27028

//Estas funciones dejan los valores en las variables new__ y el tamaño en las variables new__Size
int check_lon_prec(double *lon, int ncIdPrec, int precId,double *newLon, int *newLonSize);
int check_lon_tg(double *lon, int ncIdTg, int tgId,double *newLon, int *newLonSize);
int check_lat_prec(double *lat, int ncIdPrec, int precId,double *newLat, int *newLatSize);
int check_lat_tg(double *lat, int ncIdTg, int tgId,double *newLat, int *newLatSize);

//Estas funciones comprueban, interpolan y escriben en el nuevo fichero.
int check_and_write_prec(double *oldLon, double *oldLat,int ncIdPrec, int precId, int newNcId, int newPrecId, int newLatSize, int newLonSize);
int check_and_write_tg(double *oldLon, double *oldLat,int ncIdTg, int tgId, int newNcId, int newTgId, int newLatSize, int newLonSize);
 /*
  El formato del fichero sin optimizar es TIME, LATITUD y LONGITUD, en ese orden de acceso
  Para que ambos dataset tengan las mismas coordenadas, se ha decidido solo quedarse con las 
coordenadas de la precicipitación, que tiene menos coordenadas tras el filtrado. Para eliminar 
esta opción descomentar la parte de la temperatura
*/
int main(){
  int i,ncIdPrec, ncIdTg, lonId, latId, timeId, newLonSize=0, newLatSize=0, precId, tgId, retval;
  int newNcId, newLonDim, newLatDim, newTimeDim, dimIds[3], newPrecId, newTgId, newLonId, newLatId, newTimeId; //Ida for new file
  double *newLon, *newLat;
  double lat[LAT];
  double lon[LON];
  double time[TIME];
  struct timeval t, t2;
  double segundos;
  gettimeofday(&t, NULL);

  //Primero vamos a optimizar la precicpitación
  //Lectura de datos
  printf("**Creando el fichero limpio de precipitacion**\n");
  if((retval=nc_create("copernicus_prec_data_clean.nc", NC_CLOBBER, &newNcId)))
    ERR(retval);
  printf("**Fichero creado con éxito**\n");
  printf("**Lectura de datos**\n");
  if((retval=nc_open(FILE_NAME_PREC, NC_NOWRITE, &ncIdPrec)))
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
  if((retval=nc_inq_varid(ncIdPrec, "rr", &precId)))
    ERR(retval);
  printf("**Comprobando la Longitud**\n");
  //Reservamos memoria para la longitud necesaria
  if((newLon = calloc(sizeof(double),LON)) == NULL)
    MEMERR("Reserva de newLon");
  check_lon_prec(lon,ncIdPrec,precId, newLon, &newLonSize);
  printf("**Longitud Comprobada con éxito**\n");
  printf("**Declarando Longitud y tiempo en el nuevo fichero**\n");
  if((retval=nc_def_dim(newNcId, "time", TIME, &newTimeDim)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "time", NC_DOUBLE, 1, &newTimeDim, &newTimeId)))
    ERR(retval);
  if((retval=nc_def_dim(newNcId, "longitude",newLonSize, &newLonDim)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "longitude", NC_DOUBLE, 1, &newLonDim, &newLonId)))
    ERR(retval);
  if((newLat = calloc(sizeof(double), LAT)) == NULL)
    MEMERR("Reserva de newLon");
  printf("**Comprobando Latitud**\n");
  check_lat_prec(lat,ncIdPrec,precId, newLat, &newLatSize);
  printf("**Latitud Comprobada con éxito**\n");
  printf("**Declarando Latitud  y prec en el nuevo fichero, y escribiendo Lat,lon y time**\n");
  if((retval=nc_def_dim(newNcId, "latitude", newLatSize, &newLatDim)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "latitude", NC_DOUBLE, 1, &newLatDim, &newLatId)))
    ERR(retval);
  //Declaramos la variable precipitacion ahora que tenemos latitud, longitud y tiempo
  dimIds[0]=newLatDim;
  dimIds[1]=newLonDim;
  dimIds[2]=newTimeDim;
  if((retval=nc_def_var(newNcId, "rr", NC_FLOAT, 3, dimIds, &newPrecId)))
    ERR(retval);
  if((retval = nc_enddef(newNcId)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newTimeId, time)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLonId, newLon)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLatId, newLat)))
    ERR(retval);
  //NO Liberamos la latitud y longitud nuevas para ponerlas en el fichero de temperatura

  printf("**Escritura Correcta**\n");
  check_and_write_prec(lon, lat, ncIdPrec, precId, newNcId, newPrecId, newLatSize, newLonSize);
  if((retval = nc_close(newNcId)))
    ERR(retval);
  printf("**Escritura de precipitacion en fichero nuevo completado**\n");

  //Ahora la temperatura
  //Lectura de datos
  printf("**Creando el fichero limpio de temperatura**\n");
  if((retval=nc_create("copernicus_tg_data_clean.nc", NC_CLOBBER, &newNcId)))
    ERR(retval);
  printf("**Fichero creado con éxito**\n");
  printf("**Lectura de datos**\n");
  if((retval=nc_open(FILE_NAME_TG, NC_NOWRITE, &ncIdTg)))
    ERR(retval);
  /* Comentado para que los dataset tengan las mismas coordenadas
  if((retval=nc_inq_varid(ncIdTg, "longitude", &lonId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdTg, lonId, lon)))
    ERR(retval);
  if((retval=  nc_inq_varid(ncIdTg, "latitude", &latId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdTg, latId, lat)))
    ERR(retval);
  if((retval= nc_inq_varid(ncIdTg, "time", &timeId)))
    ERR(retval);
  if((retval=nc_get_var_double(ncIdTg, timeId, time)))
    ERR(retval);
  */
  if((retval=nc_inq_varid(ncIdTg, "tg", &tgId)))
    ERR(retval);
  /*
  printf("**Comprobando la Longitud**\n");
  //Reservamos memoria para la longitud necesaria
  if((newLon = calloc(sizeof(double),LON)) == NULL)
    MEMERR("Reserva de newLon");
  check_lon_tg(lon,ncIdTg,tgId, newLon, &newLonSize);
  printf("**Longitud Comprobada con éxito**\n");
  */
  printf("**Declarando Longitud y tiempo en el nuevo fichero**\n");
  if((retval=nc_def_dim(newNcId, "time", TIME, &newTimeDim)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "time", NC_DOUBLE, 1, &newTimeDim, &newTimeId)))
    ERR(retval);
  if((retval=nc_def_dim(newNcId, "longitude",newLonSize, &newLonDim)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "longitude", NC_DOUBLE, 1, &newLonDim, &newLonId)))
    ERR(retval);
  /*
  if((newLat = calloc(sizeof(double), LAT)) == NULL)
    MEMERR("Reserva de newLon");
  printf("**Comprobando Latitud**\n");
  check_lat_tg(lat,ncIdTg,tgId, newLat, &newLatSize);
  printf("**Latitud Comprobada con éxito**\n");
  */
  printf("**Declarando Latitud  y tg en el nuevo fichero, y escribiendo Lat,lon y time**\n");
  if((retval=nc_def_dim(newNcId, "latitude", newLatSize, &newLatDim)))
    ERR(retval);
  if((retval=nc_def_var(newNcId, "latitude", NC_DOUBLE, 1, &newLatDim, &newLatId)))
    ERR(retval);
  //Declaramos la variable temperatura ahora que tenemos latitud, longitud y tiempo
  dimIds[0]=newLatDim;
  dimIds[1]=newLonDim;
  dimIds[2]=newTimeDim;
  if((retval=nc_def_var(newNcId, "tg", NC_SHORT, 3, dimIds, &newTgId)))
    ERR(retval);
  if((retval = nc_enddef(newNcId)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newTimeId, time)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLonId, newLon)))
    ERR(retval);
  if ((retval = nc_put_var_double(newNcId, newLatId, newLat)))
    ERR(retval);
  //Liberamos la latitud y longitud nuevas
  free(newLat);
  free(newLon);
  printf("**Escritura Correcta**\n");
  check_and_write_tg(lon, lat, ncIdTg, tgId, newNcId, newTgId, newLatSize, newLonSize);
  if((retval = nc_close(newNcId)))
    ERR(retval);
  printf("**Escritura de Temperatura en fichero nuevo completado**\n");

  printf("**SUCCESS!!**\n");
  gettimeofday(&t2, NULL);
  segundos = (((t2.tv_usec - t.tv_usec)/1000000.0f)  + (t2.tv_sec - t.tv_sec));
 
  printf("\n *******> Duracion total: %f segundos\n\n", segundos);
  nc_close(ncIdTg);
  return 0;
}

int check_lon_prec(double *lon, int ncIdPrec, int precId,double *newLon, int *newLonSize){
  int i,j,k,retval, thread_id;
  double *partPrecVal;
  size_t *startp, *countp;
  unsigned long limite =(LAT*TIME) - ((LAT*TIME/100)*5);
  unsigned long *notValid = calloc(sizeof(unsigned long), LON);
  *(newLonSize)=0;
  //Reservamos memoria para todos los dias de una latitud
  if((partPrecVal = calloc(sizeof(double), LAT)) == NULL)
    MEMERR("Reserva de precVal para LAT");
  startp=(size_t *)malloc(NDIMS*sizeof(size_t));
  countp=(size_t *)malloc(NDIMS*sizeof(size_t));
  *countp=1; //Leemos todos los dias
  *(countp+1)=LAT; //todos los datos de la latitud
  *(countp+2)=1; //un dato para la longitud que queremos comprobar
  *(startp+1)=0;//Siempre empezamos a leer desde la primera latitud
    
  //Por cada LON, para cada día nos traemos los datos de cada LAT, aunque no son contiguos porque vamos a leer por LAT.
  for(i=0;i<TIME;i++){
    if(i == TIME/4)
      printf("**25%% del calculo completado**\n");
    else if(i == TIME/2)
      printf("**50%% del calculo completado**\n");
    else if(i == ((TIME/4)*3))
      printf("**75%% del calculo completado**\n");
    *startp = i;
    for(j=0;j<LON;j++){
      *(startp+2)= j;
      
      if((retval=nc_get_vara_double(ncIdPrec,precId,startp,countp, partPrecVal))){
	printf("error get prec while checking lon\n");
	ERR(retval);
      }
      for(k=0; k<LAT;k++){
	 if(*(partPrecVal + k) < 0.0 || *(partPrecVal + k) > 300.0){
	  *(notValid+j) = *(notValid+j)+1;
	 }
      }
    }
  }
  for(i=0;i<LON;i++){
    if(*(notValid+i) < limite){
      *(newLon+*(newLonSize)) = *(lon+i);
      *(newLonSize) = *(newLonSize) +1;
    } else {
      *(lon+i) = -9999.0;
    }
  }
  free(notValid);
  free(partPrecVal);
  free(startp);
  free(countp);
  return 0;
}
int check_lon_tg(double *lon, int ncIdTg, int tgId,double *newLon, int *newLonSize){
  int i,j,k,retval;
  short *partTgVal;
  size_t *startp, *countp;
  unsigned long limite =(LAT*TIME) - ((LAT*TIME/100)*5);
  unsigned long *notValid = calloc(sizeof(unsigned long), LON);
  *(newLonSize)=0;
  //Solo reservamos la memoria necesaria para comprobar un solo día de todas las latitudes
  if((partTgVal = calloc(sizeof(short), LAT)) == NULL)
    MEMERR("Reserva de tgVal para LAT");
   
  startp=(size_t *)malloc(NDIMS*sizeof(size_t));
  countp=(size_t *)malloc(NDIMS*sizeof(size_t));
  *countp=1; //Leemos un dato de del dia
  *(countp+1)=LAT; //todos los datos de la latitud
  *(countp+2)=1; //un dato para la longitud que queremos comprobar
  *(startp+1)=0;//Siempre empezamos a leer desde la primera latitud

  //Por cada LON, para cada día nos traemos los datos de cada LAT, aunque no son contiguos porque vamos a leer por LAT.
  for(i=0;i<TIME;i++){
    if(i == TIME/4)
      printf("**25%% del calculo completado**\n");
    else if(i == TIME/2)
      printf("**50%% del calculo completado**\n");
    else if(i == ((TIME/4)*3))
      printf("**75%% del calculo completado**\n");
    *startp = i;
    for(j=0;j<LON;j++){
      *(startp+2)= j;
      if((retval=nc_get_vara_short(ncIdTg,tgId,startp,countp, partTgVal))){
	printf("error get tg while checking lon\n");
	ERR(retval);
      }
      for(k=0; k<LAT;k++){
	if(*(partTgVal + k) < -8500 || *(partTgVal + k) > 8500)
	  *(notValid+j) = *(notValid+j)+1;
      }
    } 
  }
  for(i=0;i<LON;i++){
    if(*(notValid+i) < limite){
      *(newLon+*(newLonSize)) = *(lon+i);
      *(newLonSize) = *(newLonSize) +1;
    } else {
      *(lon+i) = -9999.0;
    }
  }
  free(notValid);
  free(partTgVal);
  free(startp);
  free(countp);
  return 0;
}
int check_lat_prec(double *lat, int ncIdPrec, int precId,double *newLat, int *newLatSize){
  int i,j,k,retval, thread_id;
  double *partPrecVal;
  size_t *startp, *countp;
  unsigned long limite =(LON*TIME) - ((LON*TIME/100)*5);
  unsigned long *notValid = calloc(sizeof(unsigned long), LAT);
  
  *(newLatSize)=0;
  //Solo reservamos la memoria necesaria para comprobar un solo día de todas las latitudes
  if((partPrecVal = (double *) calloc(sizeof(double), LON)) == NULL)
    MEMERR("Reserva de precVal para LAT");
   
  startp=(size_t *)malloc(NDIMS*sizeof(size_t));
  countp=(size_t *)malloc(NDIMS*sizeof(size_t));
  *countp=1; //Leemos un dato de del dia
  *(countp+1)=1; //un dato latitud
  *(countp+2)=LON; //todos los datos de LON
  *(startp+2)=0;//Siempre empezamos a leer desde la primera longitud
  //Por cada LAT, para cada día nos traemos los datos de cada LON, por tanto los bloques son contiguos
  for(i=0;i<TIME;i++){
    if(i == TIME/4)
      printf("**25%% del calculo completado**\n");
    else if(i == TIME/2)
      printf("**50%% del calculo completado**\n");
    else if(i == ((TIME/4)*3))
      printf("**75%% del calculo completado**\n");
    *(startp)= i;
    for(j=0;j<LAT;j++){
      *(startp+1) = j;
      if((retval=nc_get_vara_double(ncIdPrec,precId,startp,countp, partPrecVal))){
	printf("error get prec while checking lon\n");
	ERR(retval);
      }
      for(k=0; k<LON;k++){
	if(*(partPrecVal + k) < 0.0 || *(partPrecVal + k) > 300.0)
	  *(notValid+j) = *(notValid+j)+1;
      }
    }
  }
  for(i=0;i<LAT;i++){
    if(*(notValid+i) < limite){
      *(newLat+*(newLatSize)) = *(lat+i);
      *(newLatSize) = *(newLatSize) +1;
    } else {
      *(lat+i) = -9999.0;
    }
  }
  free(notValid);
  free(partPrecVal);
  free(startp);
  free(countp);
  return 0;
}
int check_lat_tg(double *lat, int ncIdTg, int tgId,double *newLat, int *newLatSize){
  int i,j,k,retval;
  short *partTgVal;
  size_t *startp, *countp;
  unsigned long limite =(LON*TIME) - ((LON*TIME/100)*5);
  unsigned long *notValid=calloc(sizeof(unsigned long), LAT);
  *(newLatSize)=0;
  //Solo reservamos la memoria necesaria para comprobar un solo día de todas las latitudes
  if((partTgVal = calloc(sizeof(short), LON)) == NULL)
    MEMERR("Reserva de tgVal para LAT");
   
  startp=(size_t *)malloc(NDIMS*sizeof(size_t));
  countp=(size_t *)malloc(NDIMS*sizeof(size_t));
  *countp=1; //Leemos un dato de del dia
  *(countp+1)=1; //un dato latitud
  *(countp+2)=LON; //todos los datos de LON
  *(startp+2)=0;//Siempre empezamos a leer desde la primera longitud

  //Por cada LAT, para cada día nos traemos los datos de cada LON, por tanto los bloques son contiguos.
  for(i=0;i<TIME;i++){
    if(i == TIME/4)
      printf("**25%% del calculo completado**\n");
    else if(i == TIME/2)
      printf("**50%% del calculo completado**\n");
    else if(i == ((TIME/4)*3))
      printf("**75%% del calculo completado**\n");
    *startp = i;
    for(j=0;j<LAT;j++){
      *(startp+1)= j;
      if((retval=nc_get_vara_short(ncIdTg,tgId,startp,countp, partTgVal))){
	printf("error get tg while checking lon\n");
	ERR(retval);
      }
      for(k=0; k<LON;k++){
	if(*(partTgVal + k) < -8500 || *(partTgVal + k) > 8500)
	  *(notValid+j) = *(notValid+j)+1;
      }
    }
  }
  for(i=0;i<LAT;i++){
    if(*(notValid+i) < limite){
      *(newLat+(*newLatSize)) = *(lat+i);
      *(newLatSize) = *(newLatSize) +1;
    } else {
      *(lat+i) = -9999.0;
    }
  }
  free(notValid);
  free(partTgVal);
  free(startp);
  free(countp);
  return 0;
}
/* Debido a que se puede llegar a leer un gran volumen de datos, y vamos a transformar los datos de forma que esten de forma optimizada
para acceder de forma temporal, se ha decidido dividir las escrituras. Como las escrituras deben ser directas, se ha establecido una variable
global que establece en cuantas veces se dividirán las escrituras. Divide el número de latitudes que se escribirán a la vez. NUM_ESCR_PREC
*/
int check_and_write_prec(double *oldLon, double *oldLat,int ncIdPrec, int precId, int newNcId, int newPrecId, int newLatSize, int newLonSize){
  int i,j,k,z, utilLatToWrite=0, writeLatIndex=0;
  int retval;
  unsigned long newPrecIndex =0, despl_total;
  size_t *startp,*countp, *writeStartp, *writeCountp;
  double *partPrecVal;
  float *escrPrecVal;
  int lastLatThisIter=0;
  int numLatPerIter= LAT/NUM_ESCR_PREC;
  //Vamos a leer 3 Lat a la vez para poder interpolar con los datos de alrededor.
  if((partPrecVal=malloc(sizeof(double)*LON*3)) == NULL)
    MEMERR("Reserva de tgVal para LAT");
  //Reservamos necesario para el máximo número posible ha escribir
  if((escrPrecVal=malloc(sizeof(float)*((numLatPerIter+LAT%NUM_ESCR_PREC)*newLonSize*TIME))) == NULL)
    MEMERR("Reserva de escrPrecVal para LAT");
  
  //Para lecturas. Orden del fichero de lectura TIME -> LAT -> LON
  startp=(size_t *)malloc(NDIMS*sizeof(size_t));
  countp=(size_t *)malloc(NDIMS*sizeof(size_t));
  *countp=1; //Leemos un dato de del dia
  *(countp+1)=3; // para poder interpolar con los datos a la derecha e izquierda del punto
  *(countp+2)=LON; //todos los datos de LON
  *(startp+2)=0;//Siempre empezamos a leer desde la primera longitud
  //Para escrituras. Orden del fichero nuevo LAT -> LON -> TIME
  writeStartp=(size_t *)malloc(NDIMS*sizeof(size_t));
  writeCountp=(size_t *)malloc(NDIMS*sizeof(size_t));
  *(writeCountp+2)=TIME; // Escribimos todos los datos de un dia
  *(writeCountp+1)=newLonSize; //Escribimos todos los datos de LON
  *(writeStartp+1)=0;//Siempre empezamos a escribir desde la primera LON
  *(writeStartp+2)=0;//Siempre empezamos a escribir desde el primer dia
  printf("El fichero se ha dividido en %d para escribir\n",NUM_ESCR_PREC);
  for(z=0; z<NUM_ESCR_PREC; z++){
  
    utilLatToWrite=0;
    lastLatThisIter=numLatPerIter*(z+1);//La ultima LAT de esta iter
    if(z==NUM_ESCR_PREC-1){//Si es la ultima hacemos las sobrantes
      lastLatThisIter=LAT;
    }
    for(i=0; i<TIME; i++){
      *startp=i;
      utilLatToWrite=0;
      //Marcamos las latitudes que nos tocan esta iter
      for(j=z*numLatPerIter; j<lastLatThisIter; j++){
	if(*(oldLat+j) == -9999.0)
	  continue;
	utilLatToWrite++;
	*(startp+1)=j-1;
	//Leemos la precipitación. Vamos a mirar la longitud del medio, por lo que nos traemos los de alrededor.
	if((retval=nc_get_vara_double(ncIdPrec,precId,startp,countp, partPrecVal))){
	  printf("error get prec while checking lon\n");
	  ERR(retval);
	}
	//La latitud extrema izquierda
	if(j=1){
	  newPrecIndex=i;//Iniciamos el indice por el dia por el que vamos
	  for(k=0; k<LON; k++){
	    if(*(oldLon+k) == -9999.0)
	      continue;
	    //Interpolamos y copiamos los datos útiles
	    despl_total = k; //Nos movemos a la latitud del medio
	    if(*(partPrecVal+despl_total) < 0){ //En caso de que sea nulo, copiamos el valor de un lugar cercano
	      if(k>0 && *(partPrecVal+despl_total-1) >= 0.0) //Arriba
		*(partPrecVal+despl_total) = *(partPrecVal+despl_total-1);
	      else if(k<LON-1 && *(partPrecVal+despl_total+1) >= 0.0) //Abajo
		*(partPrecVal+despl_total) = *(partPrecVal+despl_total+1);
	      else if(*(partPrecVal+despl_total+LON) >= 0.0) //Derecha
		*(partPrecVal+despl_total) = *(partPrecVal+despl_total+LON);
	    }
	    //Desplazamos el puntero para que todas las posiciones tengan sus días contiguos.
	    *(escrPrecVal+newPrecIndex)= (float) *(partPrecVal+despl_total);
	    newPrecIndex = newPrecIndex + TIME;
	  }	
	}
	//Iniciamos el indice por el dia en el que vamos(i), y en el espacio de la latitud por la que vamos(j-z*numLatIter), y cada latitud contiene newLonSize LON con tantos dias como TIME
	newPrecIndex=i+utilLatToWrite*newLonSize*TIME;
	for(k=0; k<LON; k++){
	  if(*(oldLon+k) == -9999.0)
	    continue;
	  //Interpolamos y copiamos los datos útiles
	  despl_total = k+LON; //Nos movemos a la latitud del medio
	  if(*(partPrecVal+despl_total) < 0){ //En caso de que sea nulo, copiamos el valor de un lugar cercano
	    if(k>0 && *(partPrecVal+despl_total-1) >= 0.0) //Arriba
	      *(partPrecVal+despl_total) = *(partPrecVal+despl_total-1);
	    else if(k<LON-1 && *(partPrecVal+despl_total+1) >= 0.0) //Abajo
	      *(partPrecVal+despl_total) = *(partPrecVal+despl_total+1);
	    else if(*(partPrecVal+despl_total+LON) >= 0.0) //Derecha
	      *(partPrecVal+despl_total) = *(partPrecVal+despl_total+LON);
	    else if(*(partPrecVal+despl_total-LON) >= 0.0) //Izquierda
	      *(partPrecVal+despl_total) = *(partPrecVal+despl_total-LON);
	  }
	  *(escrPrecVal+newPrecIndex)= (float) *(partPrecVal+despl_total);
	  //Desplazamos el puntero para que todas las posiciones tengan sus días contiguos.
	  newPrecIndex = newPrecIndex + TIME;
	}
	//La latidud extrema derecha
	if(j=LAT-2){
	  //Iniciamos el indice por el dia en el que vamos(i), y en el espacio de la latitud por la que vamos(j-z*numLatIter), y cada latitud contiene newLonSize LON con tantos dias como TIME
	  newPrecIndex=i+utilLatToWrite*newLonSize*TIME;
	  for(k=0; k<LON; k++){
	    if(*(oldLon+k) == -9999.0)
	      continue;
	    //Interpolamos y copiamos los datos útiles
	    despl_total = k; //Nos movemos a la latitud del medio
	    if(*(partPrecVal+despl_total) < 0){ //En caso de que sea nulo, copiamos el valor de un lugar cercano
	      if(k>0 && *(partPrecVal+despl_total-1) >= 0.0) //Arriba
		*(partPrecVal+despl_total) = *(partPrecVal+despl_total-1);
	      else if(k<LON-1 && *(partPrecVal+despl_total+1) >= 0.0) //Abajo
		*(partPrecVal+despl_total) = *(partPrecVal+despl_total+1);
	      else if(*(partPrecVal+despl_total-LON) >= 0.0) //Izquierda
		*(partPrecVal+despl_total) = *(partPrecVal+despl_total-LON);
	    }
	   //Desplazamos el puntero para que todas las posiciones tengan sus días contiguos.
	    *(escrPrecVal+newPrecIndex)= (float) *(partPrecVal+despl_total);
	     newPrecIndex = newPrecIndex+ TIME;
	  }
	
	}
      }
    }
    //El indice por el que tenemos que empezar a insertar
    *writeStartp=writeLatIndex;
    //El numero de lats que tenemos que escribir
    *writeCountp=utilLatToWrite;
    //Escribimos en el fichero
    nc_put_vara_float(newNcId, newPrecId,writeStartp, writeCountp, escrPrecVal);
    //Actualizamos el indice por el que empezaremos a escribir la proxima iter
    writeLatIndex=utilLatToWrite;
    printf("%d de %d partes del fichero escritas\n",z+1, NUM_ESCR_PREC);
  }
  free(partPrecVal);
  free(escrPrecVal);
  free(startp);
  free(countp);
  free(writeStartp);
  free(writeCountp);
  return 0;
}
/*
 Variable para dividir las escrituras. NUM_ESCR_TG
 */
  int check_and_write_tg(double *oldLon, double *oldLat,int ncIdTg, int tgId, int newNcId, int newTgId, int newLatSize, int newLonSize){
 int i,j,k,z, utilLatToWrite=0, writeLatIndex=0;
  int retval;
  unsigned long newTgIndex =0, despl_total;
  size_t *startp,*countp, *writeStartp, *writeCountp;
  short *partTgVal;
  short *escrTgVal;
  int lastLatThisIter=0;
  int numLatPerIter= LAT/NUM_ESCR_TG;
  //Vamos a leer 3 Lat a la vez para poder interpolar con los datos de alrededor.
  if((partTgVal=malloc(sizeof(short)*LON*3)) == NULL)
    MEMERR("Reserva de tgVal para LAT");
  //Reservamos necesario para el máximo número posible ha escribir
  if((escrTgVal=malloc(sizeof(short)*((numLatPerIter+LAT%NUM_ESCR_TG)*newLonSize*TIME))) == NULL)
    MEMERR("Reserva de escrTgVal para LAT");
  
  //Para lecturas. Orden del fichero de lectura TIME -> LAT -> LON
  startp=(size_t *)malloc(NDIMS*sizeof(size_t));
  countp=(size_t *)malloc(NDIMS*sizeof(size_t));
  *countp=1; //Leemos un dato de del dia
  *(countp+1)=3; // para poder interpolar con los datos a la derecha e izquierda del punto
  *(countp+2)=LON; //todos los datos de LON
  *(startp+2)=0;//Siempre empezamos a leer desde la primera longitud
  //Para escrituras. Orden del fichero nuevo LAT -> LON -> TIME
  writeStartp=(size_t *)malloc(NDIMS*sizeof(size_t));
  writeCountp=(size_t *)malloc(NDIMS*sizeof(size_t));
  *(writeCountp+2)=TIME; // Escribimos todos los datos de un dia
  *(writeCountp+1)=newLonSize; //Escribimos todos los datos de LON
  *(writeStartp+1)=0;//Siempre empezamos a escribir desde la primera LON
  *(writeStartp+2)=0;//Siempre empezamos a escribir desde el primer dia
  printf("El fichero se ha dividido en %d para escribir\n",NUM_ESCR_TG);
  for(z=0; z<NUM_ESCR_TG; z++){
    utilLatToWrite=0;
    lastLatThisIter=numLatPerIter*(z+1);//La ultima LAT de esta iter
    if(z==NUM_ESCR_TG-1){//Si es la ultima hacemos las sobrantes
      lastLatThisIter=LAT;
    }
    for(i=0; i<TIME; i++){
      *startp=i;
      utilLatToWrite=0;
      //Marcamos las latitudes que nos tocan esta iter
      for(j=z*numLatPerIter; j<lastLatThisIter; j++){
	if(*(oldLat+j) == -9999.0)
	  continue;
	utilLatToWrite++;
	*(startp+1)=j-1;
	//Leemos la temperatura. Vamos a mirar la longitud del medio, por lo que nos traemos los de alrededor.
	if((retval=nc_get_vara_short(ncIdTg,tgId,startp,countp, partTgVal))){
	  printf("error get tg while checking lon\n");
	  ERR(retval);
	}
	//La latitud extrema izquierda
	if(j=1){
	  newTgIndex=i;//Iniciamos el indice por el dia por el que vamos
	  for(k=0; k<LON; k++){
	    if(*(oldLon+k) == -9999.0)
	      continue;
	    //Interpolamos y copiamos los datos útiles
	    despl_total = k; //Nos movemos a la latitud del medio
	    if(*(partTgVal+despl_total) < 0){ //En caso de que sea nulo, copiamos el valor de un lugar cercano
	      if(k>0 && *(partTgVal+despl_total-1) >= 0.0) //Arriba
		*(partTgVal+despl_total) = *(partTgVal+despl_total-1);
	      else if(k<LON-1 && *(partTgVal+despl_total+1) >= 0.0) //Abajo
		*(partTgVal+despl_total) = *(partTgVal+despl_total+1);
	      else if(*(partTgVal+despl_total+LON) >= 0.0) //Derecha
		*(partTgVal+despl_total) = *(partTgVal+despl_total+LON);
	    }
	    //Desplazamos el puntero para que todas las posiciones tengan sus días contiguos.
	    *(escrTgVal+newTgIndex)=  *(partTgVal+despl_total);
	    newTgIndex = newTgIndex + TIME;
	  }	
	}
	//Iniciamos el indice por el dia en el que vamos(i), y en el espacio de la latitud por la que vamos(j-z*numLatIter), y cada latitud contiene newLonSize LON con tantos dias como TIME
	newTgIndex=i+utilLatToWrite*newLonSize*TIME;
	for(k=0; k<LON; k++){
	  if(*(oldLon+k) == -9999.0)
	    continue;
	  //Interpolamos y copiamos los datos útiles
	  despl_total = k+LON; //Nos movemos a la latitud del medio
	  if(*(partTgVal+despl_total) < 0){ //En caso de que sea nulo, copiamos el valor de un lugar cercano
	    if(k>0 && *(partTgVal+despl_total-1) >= 0.0) //Arriba
	      *(partTgVal+despl_total) = *(partTgVal+despl_total-1);
	    else if(k<LON-1 && *(partTgVal+despl_total+1) >= 0.0) //Abajo
	      *(partTgVal+despl_total) = *(partTgVal+despl_total+1);
	    else if(*(partTgVal+despl_total+LON) >= 0.0) //Derecha
	      *(partTgVal+despl_total) = *(partTgVal+despl_total+LON);
	    else if(*(partTgVal+despl_total-LON) >= 0.0) //Izquierda
	      *(partTgVal+despl_total) = *(partTgVal+despl_total-LON);
	  }
	  *(escrTgVal+newTgIndex)= *(partTgVal+despl_total);
	  //Desplazamos el puntero para que todas las posiciones tengan sus días contiguos.
	  newTgIndex = newTgIndex + TIME;
	}
	//La latidud extrema derecha
	if(j=LAT-2){
	  //Iniciamos el indice por el dia en el que vamos(i), y en el espacio de la latitud por la que vamos(j-z*numLatIter), y cada latitud contiene newLonSize LON con tantos dias como TIME
	  newTgIndex=i+utilLatToWrite*newLonSize*TIME;
	  for(k=0; k<LON; k++){
	    if(*(oldLon+k) == -9999.0)
	      continue;
	    //Interpolamos y copiamos los datos útiles
	    despl_total = k; //Nos movemos a la latitud del medio
	    if(*(partTgVal+despl_total) < 0){ //En caso de que sea nulo, copiamos el valor de un lugar cercano
	      if(k>0 && *(partTgVal+despl_total-1) >= 0.0) //Arriba
		*(partTgVal+despl_total) = *(partTgVal+despl_total-1);
	      else if(k<LON-1 && *(partTgVal+despl_total+1) >= 0.0) //Abajo
		*(partTgVal+despl_total) = *(partTgVal+despl_total+1);
	      else if(*(partTgVal+despl_total-LON) >= 0.0) //Izquierda
		*(partTgVal+despl_total) = *(partTgVal+despl_total-LON);
	    }
	   //Desplazamos el puntero para que todas las posiciones tengan sus días contiguos.
	    *(escrTgVal+newTgIndex)= *(partTgVal+despl_total);
	     newTgIndex = newTgIndex+ TIME;
	  }
	
	}
      }
    }
    //El indice por el que tenemos que empezar a insertar
    *writeStartp=writeLatIndex;
    //El numero de lats que tenemos que escribir
    *writeCountp=utilLatToWrite;
    //Escribimos en el fichero
    nc_put_vara_short(newNcId, newTgId,writeStartp, writeCountp, escrTgVal);
    //Actualizamos el indice por el que empezaremos a escribir la proxima iter
    writeLatIndex=utilLatToWrite;
    printf("%d de %d partes del fichero escritas\n",z+1, NUM_ESCR_TG);
  }
  free(partTgVal);
  free(escrTgVal);
  free(startp);
  free(countp);
  free(writeStartp);
  free(writeCountp);
  return 0;
}
