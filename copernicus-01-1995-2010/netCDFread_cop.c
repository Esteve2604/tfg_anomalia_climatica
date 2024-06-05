#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define FILE_NAME_PREC "./rr_ens_mean_0.1deg_reg_1995-2010_v29.0e.nc"
#define FILE_NAME_TG "./tg_ens_mean_0.1deg_reg_1995-2010_v29.0e.nc"
#define NDIMS 3
#define LAT 465
#define LON 705
#define TIME 5844

#define SHORT 1
#define DOUBLE 0
int data_analisis(void *, double[], double[], double[], int);
int main(){
  int i,j,k, ncidp;
  size_t ldimsp, lvarsp, lattr, *indexp;
  nc_type varType, attrType;
  float varValue;
  short *tg_val;
  double *prec_val, *last_drought_code;
  int prec_varid, lat_varid, lon_varid, time_varid;
  double lat[LAT];
  double lon[LON];
  double time[TIME];
  struct timeval t, t2;
  double segundos;
  gettimeofday(&t, NULL);
  prec_val = calloc(sizeof(double), LON * LAT * TIME);
  nc_open(FILE_NAME_PREC, NC_NOWRITE, &ncidp);
  nc_inq_varid(ncidp, "longitude", &lon_varid);
  nc_get_var_double(ncidp, lon_varid, lon);
  nc_inq_varid(ncidp, "latitude", &lat_varid);
  nc_get_var_double(ncidp, lat_varid, lat);
  nc_inq_varid(ncidp, "time", &time_varid);
  nc_get_var_double(ncidp, time_varid, time);
  nc_inq_varid(ncidp, "rr", &prec_varid);
  nc_get_var_double(ncidp, prec_varid, prec_val);
  data_analisis(prec_val, lat, lon, time, DOUBLE);
  nc_close(ncidp);
  free(prec_val);
  tg_val = calloc(sizeof(short), LON * LAT * TIME);
  nc_open(FILE_NAME_TG, NC_NOWRITE, &ncidp);
  nc_inq_varid(ncidp, "longitude", &lon_varid);
  nc_get_var_double(ncidp, lon_varid, lon);
  nc_inq_varid(ncidp, "latitude", &lat_varid);
  nc_get_var_double(ncidp, lat_varid, lat);
  nc_inq_varid(ncidp, "time", &time_varid);
  nc_get_var_double(ncidp, time_varid, time);
  nc_inq_varid(ncidp, "tg", &prec_varid);
  nc_get_var_short(ncidp, prec_varid, tg_val);
  data_analisis(tg_val, lat, lon, time, SHORT);
  nc_close(ncidp);
  free(tg_val);
  gettimeofday(&t2, NULL);
  segundos = (((t2.tv_usec - t.tv_usec)/1000000.0f)  + (t2.tv_sec - t.tv_sec));
  
  printf("\n *******> Duracion total: %f segundos\n\n", segundos);
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
	  if(*((short *)val+ despl_total) < (short) -8500 || *((short *)val+ despl_total) > (short) 8500){
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
   printf("All Data: Max Lon: %f Max Lat: %f Min Lon: %f Min Lat: %f\n", lon[LON-1], lat[LAT-1], lon[0], lat[0]);
    printf("Valid: Max Lon: %f Max Lat: %f Min Lon: %f Min Lat: %f\n", max_lon, max_lat, min_lon, min_lat);
  printf("Total Data: %f Not Valid %ld, Percentage NV: %f%% Average_prec: %f Average_tg: %ld \n", total_data, not_valid,(not_valid/total_data)*100 , average, average_tg); 
  return 0;
}
