#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fcntl.h>  // Para manipular pipes
#include <sys/time.h> // Para usleep

// gcc ACOsecuencial.c -o ACOsecuencial -lm

#define NUM_HORMIGAS 100
#define MAX_ITERACIONES 10
#define ALPHA 1.0
#define BETA 2.0
#define EVAPORACION 0.5
#define Q 100.0

int NUM_PUNTOS = 0;
int PUNTO_INICIAL = 0;
int PUNTO_FINAL = 0;

// Matriz dinámica de feromonas
double **feromonas;

// Función para inicializar la matriz de feromonas
void inicializar_feromonas() 
{
    feromonas = (double **)malloc(NUM_PUNTOS * sizeof(double *));
    if (feromonas == NULL) 
    {
        printf("error al reservar matriz\n");
        perror("Error al reservar memoria para la matriz de feromonas");
    }
    for (int i = 0; i < NUM_PUNTOS; i++) 
    {
        feromonas[i] = (double *)malloc(NUM_PUNTOS * sizeof(double));
        if (feromonas[i] == NULL) 
        {
            printf("error al reservar matriz 2\n");
            perror("Error al reservar memoria para la matriz de feromonas");
        }
        for (int j = 0; j < NUM_PUNTOS; j++) 
        {
            feromonas[i][j] = 1.0; // Inicializar feromonas
        }
    }
}

// Función para liberar la memoria de la matriz de feromonas
void liberar_feromonas() 
{
    for (int i = 0; i < NUM_PUNTOS; i++) 
    {                
        free(feromonas[i]);
    }
    free(feromonas);
}

// Función para enviar un comando al servidor y recibir la respuesta
void comunicacion(char* pipe_wd, char* pipe_rd, const char *comando, const char *datos, char *respuesta) 
{
    char mensaje[256];
    snprintf(mensaje, sizeof(mensaje), "%s|%s", comando, datos);
    usleep(200); //Espera para enviar
    int write_fd = open(pipe_wd, O_WRONLY);
    if (write_fd == -1) 
    {
        perror("Error al abrir pipe_wd");
    }
    if (write(write_fd, mensaje, strlen(mensaje) + 1) == -1) 
    {
        printf("Error al escribir\n");
        perror("Error al escribir en el pipe");
    }
    close(write_fd);

    int read_fd = open(pipe_rd, O_RDONLY);
    if (read_fd == -1) 
    {
        perror("Error al abrir pipe_rd");
    }
    if (read(read_fd, respuesta, 256) == -1) 
    {
        perror("Error al leer del pipe");
        printf("Error al leer\n");
    }
    close(read_fd);
}


// Dadas las coordenadas del punto, me da el número de punto dentro del rectángulo que define el oceano
int obtener_numeracion(float y, float x, float x_min, float x_max, float y_min, float y_max, float resolucion) 
{
    // Calcular puntos por fila
    int puntos_por_fila = (int)floor((x_max - x_min) / resolucion) + 1;

    // Determinar índices del punto
    int indice_fila = (int)floor((y - y_min) / resolucion);    // Índice de la fila
    int indice_columna = (int)floor((x - x_min) / resolucion); // Índice de la columna

    // Calcular la numeración del punto
    int numero_punto = (indice_fila * puntos_por_fila) + indice_columna;

    return numero_punto;
}


// Dado un punto (numero), me devuelve sus coordenadas de longitud y latitud (x, y)
void obtener_coordenadas(int numero, float x_min, float x_max, float y_min, float y_max, float resolucion, float *x, float *y) 
{
    // Calcular puntos por fila
    int puntos_por_fila = (int)floor((x_max - x_min) / resolucion) + 1;

    // Calcular fila y columna del punto
    int fila = floor(numero / puntos_por_fila);
    int columna = numero % puntos_por_fila;

    // Calcular coordenadas
    *x = x_min + (columna * resolucion);
    *y = y_min + (fila * resolucion);
}

// Función para obtener el tiempo entre dos puntos
float tiempo(int origen, int destino, char* pipe_wd, char* pipe_rd) 
{
    float time = 0.0;
    char datos[100], respuesta[256] = "";
    snprintf(datos, sizeof(datos), "%d %d", origen, destino);
    comunicacion(pipe_wd, pipe_rd, "calcular_tiempo", datos, respuesta);

    // Verificar la respuesta
    if (strncmp(respuesta, "Error", 5) == 0) 
    {
        fprintf(stderr, "Error del servidor: %s\n", respuesta);
        return -1.0;
    }
    time = atof(respuesta);
    if(isnan(time))
        return -1.0;
    return time;
}

// Función para calcular la longitud de una ruta
double calcular_longitud_ruta(int *ruta, char* pipe_wd, char* pipe_rd) 
{
    double longitud = 0.0;
    int i = 0;

    while (i < NUM_PUNTOS-1 && ruta[i] != -1 && ruta[i+1] != -1) // Descarta -1
    { 
        longitud += tiempo(ruta[i], ruta[i + 1], pipe_wd, pipe_rd);
        i++;
    }
    
    return longitud;
}

// Función para seleccionar la siguiente ciudad
int seleccionar_ciudad(int ciudad_actual, int *visitado, char* pipe_wd, char* pipe_rd) 
{
    double *probabilidades = (double *)calloc(NUM_PUNTOS, sizeof(double));
    double suma = 0.0;
    int ultimo_no_visitado = -1;
    float t = 0.0;
    double heuristica = 0.0;

    for (int i = 0; i < NUM_PUNTOS; i++) 
    {
        if (!visitado[i]) 
        {
            t = tiempo(ciudad_actual, i, pipe_wd, pipe_rd);
            heuristica = 1.0 / (t + 1e-6);
            probabilidades[i] = pow(feromonas[ciudad_actual][i], ALPHA) * pow(heuristica, BETA);
            suma += probabilidades[i];
            ultimo_no_visitado = i;
        } 
        else 
            probabilidades[i] = 0.0;
    }

    // Normalizar probabilidades
    for (int i = 0; i < NUM_PUNTOS; i++) 
    {
        probabilidades[i] /= suma;
        if (i > 0) 
            probabilidades[i] += probabilidades[i - 1];
    }
    probabilidades[ultimo_no_visitado] = 1.0;

    // Seleccionar ciudad
    double r = 0.0;
    r = rand() / (double)RAND_MAX;
    for (int i = 0; i < NUM_PUNTOS; i++) 
    {
        if (r <= probabilidades[i]) 
        {
            free(probabilidades);
            return i;
        }
    }

    free(probabilidades);
    return -1;
}

// Función para actualizar las feromonas
void actualizar_feromonas(int **rutas, double longitudes[NUM_HORMIGAS]) 
{
    for (int i = 0; i < NUM_PUNTOS; i++) 
    {
        for (int j = 0; j < NUM_PUNTOS; j++) 
        {
            feromonas[i][j] *= (1 - EVAPORACION);
        }
    }
    for (int k = 0; k < NUM_HORMIGAS; k++) 
    {
        double contribucion = Q / longitudes[k];
        for (int i = 0; i < NUM_PUNTOS - 1 && rutas[k][i] != -1 && rutas[k][i + 1] != -1; i++) 
        {
            int c1 = rutas[k][i], c2 = rutas[k][i + 1];
            feromonas[c1][c2] += contribucion;
            feromonas[c2][c1] += contribucion;
        }
    }
}

// Función principal
int main() 
{

    float lat_i, lon_i, lat_f, lon_f; // Latitudes y longitudes del punto inicial y del punto final en mi trayectoria
    lat_i = -3;
    lon_i = -1;
    lat_f = -1;
    lon_f = 2;
    float lat_i_mapa, lon_i_mapa, lat_f_mapa, lon_f_mapa, resolucion_mapa; // Latitudes y longitudes que delimitan el oceano


    // Abrir pipes para comunicación
    char* pipe_wd = "pipe_wd";
    char* pipe_rd = "pipe_rd";


    srand(time(NULL));

    char respuesta[256] = "";
    comunicacion(pipe_wd, pipe_rd, "num_puntos", "", respuesta);
    NUM_PUNTOS = atoi(respuesta);
    printf("Numero de puntos del mapa: %d\n", NUM_PUNTOS);

    comunicacion(pipe_wd, pipe_rd, "limites", "", respuesta);
    sscanf(respuesta, "%f %f %f %f %f", &lat_i_mapa, &lon_i_mapa, &lat_f_mapa, &lon_f_mapa, &resolucion_mapa);
    printf("Limites del mapa recibidos del servidor: %.2f, %.2f, %.2f, %.2f, %.2f\n",
           lat_i_mapa, lon_i_mapa, lat_f_mapa, lon_f_mapa, resolucion_mapa);

    PUNTO_INICIAL = obtener_numeracion(lat_i, lon_i, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);
    PUNTO_FINAL = obtener_numeracion(lat_f, lon_f, lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa);

    printf("Punto inicial y punto final de la trayectoria: %d, %d\n", PUNTO_INICIAL, PUNTO_FINAL);

    inicializar_feromonas();

    int *mejor_ruta = (int*)malloc(NUM_PUNTOS * sizeof(int));

    double mejor_longitud = INFINITY;

    for (int iter = 0; iter < MAX_ITERACIONES; iter++) 
    {
        int **rutas = (int **)malloc(NUM_HORMIGAS * sizeof(int *));
        if (rutas == NULL) 
        {
            fprintf(stderr, "Error: fallo en malloc para rutas\n");
            exit(EXIT_FAILURE);
        }

        double longitudes[NUM_HORMIGAS];

        for (int k = 0; k < NUM_HORMIGAS; k++) 
        {
            rutas[k] = (int *)malloc(NUM_PUNTOS * sizeof(int));
            if (rutas[k] == NULL) 
            {
                fprintf(stderr, "Error: fallo en malloc para rutas[%d]\n", k);
                exit(EXIT_FAILURE);
            }

            memset(rutas[k], -1, NUM_PUNTOS * sizeof(int));

            int *visitado = (int *)calloc(NUM_PUNTOS, sizeof(int));
            rutas[k][0] = PUNTO_INICIAL;
            visitado[PUNTO_INICIAL] = 1;

            for (int i = 1; i < NUM_PUNTOS && rutas[k][i - 1] != PUNTO_FINAL; i++) 
            {
                rutas[k][i] = seleccionar_ciudad(rutas[k][i - 1], visitado, pipe_wd, pipe_rd);
                visitado[rutas[k][i]] = 1;
            }
            longitudes[k] = calcular_longitud_ruta(rutas[k], pipe_wd, pipe_rd);

            if (longitudes[k] < mejor_longitud)
            {
                mejor_longitud = longitudes[k];
                memcpy(mejor_ruta, rutas[k], NUM_PUNTOS * sizeof(int));
            }
            free(visitado);
        }

        actualizar_feromonas(rutas, longitudes);

        for (int k = 0; k < NUM_HORMIGAS; k++) 
        {
            free(rutas[k]);
        }
        free(rutas);

        printf("Iteración %d: Mejor longitud = %.2f\n", iter + 1, mejor_longitud);
    }

    printf("Mejor ruta: \n");

    float lat, lon;
    for (int i = 0; i < NUM_PUNTOS && mejor_ruta[i] != -1; i++) 
    {
        printf("%d ", mejor_ruta[i]);
        obtener_coordenadas(mejor_ruta[i], lon_i_mapa, lon_f_mapa, lat_i_mapa, lat_f_mapa, resolucion_mapa, &lon, &lat);
        printf("Latitud: %f, Longitud: %f\n", lat, lon);
    }
    
    printf("\nTiempo de la mejor ruta: %.2f\n", mejor_longitud);

    liberar_feromonas();
    free(mejor_ruta);
    comunicacion(pipe_wd, pipe_rd, "salir", "", respuesta);

    return 0;
}