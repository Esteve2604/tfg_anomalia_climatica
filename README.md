Implementación del TFG Implementación y Evaluación de un Algoritmo de Altas Prestaciones para el Cálculo de la Anomalía Climática de Esteban Aspe Ruizç
En cada una de las Carpetas, se encuentran las distintas implementaciones realizadas en el TFG.

Para ejecutarlos, es necesario tener los dataset de temperatura y precipitación de Copernicus https://cds.climate.copernicus.eu/cdsapp#!/dataset/10.24381/cds.151d3ec6?tab=overview con el periodo y la precisión indicada en el nombre de la carpeta dentro de esa carpeta.
El 'product type' a seleccionar en la web es Ensemble Mean. 
En la Carpeta 2011-2019 se utilizó una precisión de 0.25dg.

Todos los programas tienen definidos como contantes las rutas de los ficheros a utilizar o el número de hilos con los que ejecutar, por lo puede ser necesario cambiarlas.
Para compilar cada programa, es necesario tener descargada la librería de netCDF (En Ubuntu: sudo apt install libnetcdf-dev). 
Para compilar los programas con extensión .c, se pueden compilar con gcc con los flags -O3 -lm -lnetcdf -fopenmp. 
Para los programas con extensión .cu, se puede compilar con el compilador de nvidia nvcc con los flags -O3 -lm -lnetcdf -Xcompiler -fopenmp. Si las GPUs utilizadas son antiguas, debe indicarse su compute_capability con -arch. 
Las versiones asíncronas pueden dar problemas con versiones antiguas de CUDA. Se ha probado dichas versiones con CUDA de hpc_sdk versión 24.1.
