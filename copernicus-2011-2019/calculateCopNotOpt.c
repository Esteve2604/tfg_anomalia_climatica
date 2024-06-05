#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define FILE_NAME_PREC "./rr_ens_mean_0.25deg_reg_2011-2019_v21.0e.nc"
#define FILE_NAME_TG "./tg_ens_mean_0.25deg_reg_2011-2019_v21.0e.nc"
#define NDIMS 4
#define LAT 201
#define LON 464
#define TIME 3287

#define SHORT 1
#define DOUBLE 0
int data_analisis(void *, double[], double[], double[], int);
double getDayLength(int month);
double calculateDroughtCode(short temperature, double rainfall, double droughtCode, int month);

int main(){
  int i,j,k, ncIdPrec, ncIdTg;
  size_t ldimsp, lvarsp, lattr, *indexp;
  nc_type varType, attrType;
  float varValue;
  short *tg_val, *part_tg_val[LAT*LON];
  double *prec_val, *part_prec_val[1024], *last_drought_code;
  int despl_time, despl_lat, despl_total;
  int precId, tgId, lat_varid, lon_varid, time_varid;
  double lat[LAT];
  double lon[LON];
  double time[TIME];
  struct timeval t, t2;
  double segundos;
  gettimeofday(&t, NULL);
  //prec_val = calloc(sizeof(double), LON * LAT * TIME);
  last_drought_code = calloc(sizeof(double), LON * LAT);
  nc_open(FILE_NAME_PREC, NC_NOWRITE, &ncIdPrec);
  nc_open(FILE_NAME_TG, NC_NOWRITE, &ncIdTg);
  /*nc_inq_varid(ncidp, "longitude", &lon_varid);
  nc_get_var_double(ncidp, lon_varid, lon);
  nc_inq_varid(ncidp, "latitude", &lat_varid);
  nc_get_var_double(ncidp, lat_varid, lat);
  nc_inq_varid(ncidp, "time", &time_varid);
  nc_get_var_double(ncidp, time_varid, time);
  
  nc_get_var_double(ncidp, prec_varid, prec_val);
  // data_analisis(prec_val, lat, lon, time, DOUBLE);
  tg_val = calloc(sizeof(short), LON * LAT * TIME);
  
  nc_inq_varid(ncidp, "longitude", &lon_varid);
  nc_get_var_double(ncidp, lon_varid, lon);
  nc_inq_varid(ncidp, "latitude", &lat_varid);
  nc_get_var_double(ncidp, lat_varid, lat);
  nc_inq_varid(ncidp, "time", &time_varid);
  nc_get_var_double(ncidp, time_varid, time);
  
  nc_get_var_short(ncidp, prec_varid, tg_val);
  data_analisis(tg_val, lat, lon, time, SHORT);*/
  nc_inq_varid(ncIdPrec, "rr", &precId);
  nc_inq_varid(ncIdTg, "tg", &tgId);
  for(i=0;i<LAT*LON;i++){
    part_tg_val[i]=calloc(sizeof(short), TIME);
    part_prec_val[i]=calloc(sizeof(double), TIME);
  }
  /*for(i = 0; i<TIME;i++){
    despl_time = i*LAT*LON; //"Cambio de matriz" */
    for(j=0;j<LAT; j++){
      despl_lat = j * LON; //Desplazamiento Filas
      for(k=0;k<LON;k++){
	/*	despl_total = despl_time + despl_lat + k; //Desplazamiento total con desplazamiento de columna
		if(*(prec_val+ despl_total) < 0 || *(tg_val+ despl_total) < -8000 || *(tg_val+ despl_total) > 8000){
		continue;
		}*/
	
	nc_get_vars_double(ncIdPrec,precId,despl_lat+k,TIME,LAT*LON,part_prec_val[j*LON+k]);
	nc_get_vars_double(ncIdTg,tgId,despl_lat+k,TIME,LAT*LON,part_tg_val[j*LON+k]);
	calculateLastDroughtCode(part_tg_val , part_prec_val,TIME,&*(last_drought_code+despl_lat+k));
      }
    }
    //}
 gettimeofday(&t2, NULL);
 segundos = (((t2.tv_usec - t.tv_usec)/1000000.0f)  + (t2.tv_sec - t.tv_sec));

  printf("\n *******> Duracion total: %f segundos\n\n", segundos);
 free(last_drought_code);
  nc_close(ncidp);
  free(prec_val);
  nc_close(ncidp);
  free(tg_val);
  return 0;
}

int data_analisis(void *val, double lat[], double lon[], double time[], int type){
  int i, j, k;
  long not_valid=0;
  double max_lon=-9999, max_lat=-9999, min_lon=9999,min_lat=9999;
  int despl_time, despl_rlat, despl_total;
  long average_tg=0;

  float total_data=0, average=0, max_prec=0, min_prec=99999, valores_anomalos_men=0, valores_anomalos_may=0; 
  float evapotranspiration= 0.0;
  /*  for(j=0;j<LAT; j++)
    printf(" lat: %f ;",  lat[j]);
  for(k=0;k<LON;k++)
    printf(" lon:%f ; ", lon[k]);*/
  total_data=LON * LAT  * TIME;
  for(i = 0; i<TIME;i++){
    despl_time = i*LAT*LON; //"Cambio de matriz"
    for(j=0;j<LAT; j++){
      despl_rlat = j * LON; //Desplazamiento Filas
      for(k=0;k<LON;k++){
	despl_total = despl_time + despl_rlat + k; //Desplazamiento total con desplazamiento de columna
	if(type ==DOUBLE){
	  if(*((double *)val+ despl_total) < 0){
	    not_valid++;
	    continue;
	  }
	  average+=*((double *)val+despl_total);
	} else{
	  if(*((short *)val+ despl_total) < 0){
	    not_valid++;
	    continue;
	  }
	  average_tg+=*((short *)val+despl_total);
	}
	max_lon = max_lon < lon[k] ? lon[k] : max_lon;
	max_lat = max_lat < lat[j] ? lat[j] : max_lat;
	min_lon = min_lon > lon[k] ? lon[k] : min_lon;
	min_lat = min_lat > lat[j] ? lat[j] : min_lat;

      }
    }
  }
  average/=(total_data-not_valid);
  average_tg/=(total_data-not_valid);
  printf("Total Data: %f Not Valid %ld, Percentage NV: %f%% Average_prec: %f Average_tg: %ld \n", total_data, not_valid,(not_valid/total_data)*100 , average, average_tg); 
  return 0;
}
int calculateLastDroughCode(short *temperature, double *rainfall, int nDays, double *droughtCode){
  
}
double calculateDroughtCode(short temperature, double rainfall, double droughtCode, int month){
  double prevDroughtCode= droughtCode;
  double efRainfall, moisture, prevMoisture, evapotranspiration, dayLength;
  if(rainfall>2.8){
    efRainfall = 0.86 * rainfall -1.27;
    prevMoisture = 800 * exp(-droughtCode/400);
    moisture = prevMoisture + 3.937 * efRainfall;
    prevDroughtCode = 400 * log(800/moisture);
    prevDroughtCode = prevDroughtCode < 0 ? 0 : prevDroughtCode;
  }
  evapotranspiration = 0.36 * (temperature + 2.0) + getDayLength(month);
  evapotranspiration = evapotranspiration < 0 ? 0 : evapotranspiration;
  return prevDroughtCode + 0.5 * evapotranspiration;
}

double getDayLength(int day){
  int dayYear = day%365;
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
