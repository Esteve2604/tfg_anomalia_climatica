#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define FILE_NAME_PREC "./mp/copernicus_anomalia_climatica_1995_2010.nc"
#define NDIMS 3
#define LAT 437
#define LON 592
#define TIME 5844

#define SHORT 1
#define DOUBLE 0
int data_analisis(short *, double[], double[], double[]);
int main(){
  int i,j,k, ncidp;
  size_t ldimsp, lvarsp, lattr, *indexp;
  nc_type varType, attrType;
  float varValue;
  short *anomalia_val;
  int anomalia_varid, lat_varid, lon_varid, time_varid;
  double lat[LAT];
  double lon[LON];
  double time[TIME];
   struct timeval t, t2;
  double segundos;
  gettimeofday(&t, NULL);
  anomalia_val = calloc(sizeof(short), LON * LAT * TIME);
  nc_open(FILE_NAME_PREC, NC_NOWRITE, &ncidp);
  nc_inq_varid(ncidp, "longitude", &lon_varid);
  nc_get_var_double(ncidp, lon_varid, lon);
  nc_inq_varid(ncidp, "latitude", &lat_varid);
  nc_get_var_double(ncidp, lat_varid, lat);
  nc_inq_varid(ncidp, "time", &time_varid);
  nc_get_var_double(ncidp, time_varid, time);
  nc_inq_varid(ncidp, "anomaliaClimatica", &anomalia_varid);
  nc_get_var_short(ncidp, anomalia_varid, anomalia_val);
  data_analisis(anomalia_val, lat, lon, time);
  nc_close(ncidp);
  free(anomalia_val);
  gettimeofday(&t2, NULL);
  segundos = (((t2.tv_usec - t.tv_usec)/1000000.0f)  + (t2.tv_sec - t.tv_sec));
  
  printf("\n *******> Duracion total: %f segundos\n\n", segundos);

  return 0;
}

int data_analisis(short *val, double lat[], double lon[], double time[]){
  int i, j, k;
  long not_valid=0;
  double max_lon=-9999, max_lat=-9999, min_lon=9999,min_lat=9999;
  int despl_time, despl_rlat, despl_total;
  long average_anom=0;

  float total_data=0, average=0, max_prec=0, min_prec=99999, valores_anomalos_men=0, valores_anomalos_may=0; 
  float evapotranspiration= 0.0;
  total_data=LON * LAT  * TIME;
  for(i = 0; i<LAT;i++){
    despl_time = i*LON*TIME; // "Cambio de matriz"
    for(j=0;j<LON; j++){
      despl_rlat = j * TIME; //Desplazamiento Filas
      for(k=0;k<TIME;k++){
	despl_total = despl_time + despl_rlat + k; //Desplazamiento total con desplazamiento de column
	if(*(val+ despl_total) < 0 || *((short *)val+ despl_total) > (short) 99){
	  not_valid++;
	  continue;
	}
	average_anom+=*((short *)val+despl_total);
	max_lon = max_lon < lon[j] ? lon[j] : max_lon;
	max_lat = max_lat < lat[i] ? lat[i] : max_lat;
	min_lon = min_lon > lon[j] ? lon[j] : min_lon;
	min_lat = min_lat > lat[i] ? lat[i] : min_lat;
      }
    }
  }
  average_anom/=(total_data-not_valid);
   printf("All Data: Max Lon: %f Max Lat: %f Min Lon: %f Min Lat: %f\n", lon[LON-1], lat[LAT-1], lon[0], lat[0]);
    printf("Valid: Max Lon: %f Max Lat: %f Min Lon: %f Min Lat: %f\n", max_lon, max_lat, min_lon, min_lat);
  printf("Total Data: %f Not Valid %ld, Percentage NV: %f%% Average_tg: %ld \n", total_data, not_valid,(not_valid/total_data)*100 , average_anom); 
  return 0;
}
