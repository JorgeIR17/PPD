import os
import math
import time
import numpy as np
import argparse
import math
from src.data import Ocean
from src.simulation import compute_time

#python3 script11.py -3 -1 -1 2 1 1

# Variables globales para la configuración inicial
VEL_BARCO = 0  # Velocidad predeterminada del barco en m/s
RESOLUCION = 0  # Resolución predeterminada del mapa

# Inicializar la base de datos del océano
ocean_data = Ocean("data/vo.csv", "data/uo.csv")


def num_puntos(lat_i, lon_i, lat_f, lon_f):
    """
    Calcula el número de puntos en una cuadrícula definida por latitudes y longitudes.
    """
    lat_min, lat_max = min(lat_i, lat_f), max(lat_i, lat_f)
    lon_min, lon_max = min(lon_i, lon_f), max(lon_i, lon_f)
    n_puntos_lat = int((lat_max - lat_min) // RESOLUCION + 1)
    n_puntos_lon = int((lon_max - lon_min) // RESOLUCION + 1)
    return n_puntos_lat * n_puntos_lon

# Dado un punto, me devuelve sus coordenadas de latitud y longitud
def obtener_coordenadas(numero, x_min, x_max, y_min, y_max):
    # Paso 1: Calcular puntos por fila
    puntos_por_fila = int((x_max - x_min) // RESOLUCION) + 1

    # Paso 2: Calcular fila y columna del punto
    fila = numero // puntos_por_fila
    columna = numero % puntos_por_fila

    # Paso 3: Calcular coordenadas
    x = x_min + (columna * RESOLUCION)
    y = y_min + (fila * RESOLUCION)
    
    return (x, y)

def tiempo(lat_i, lon_i, lat_f, lon_f):
    """
    Calcula el tiempo necesario para viajar entre dos puntos considerando la velocidad del barco.
    """
    return compute_time(
        lat_start=lat_i,
        lon_start=lon_i,
        lat_end=lat_f,
        lon_end=lon_f,
        vel_ship=VEL_BARCO,
        ocean_data=ocean_data,
    )


def comunicacion(comando, datos, lati, loni, latf, lonf):
    """
    Procesa un comando recibido del ACO y devuelve la respuesta correspondiente.
    """
    try:
        if comando == "salir":
            return "Servidor cerrado", True
        elif comando == "limites":
            limites = f"{lati} {loni} {latf} {lonf} {RESOLUCION}"
            return (limites), False
        
        elif comando == "calcular_tiempo":
            try:
                datos = datos.rstrip('\x00')
                datos = datos.split()
                punto_inicio = float(datos[0])
                punto_fin = float(datos[1])
            except ValueError as e:
                return f"Error puntos: {e}", False

            coordenadas_inicio = obtener_coordenadas(punto_inicio, loni, lonf, lati, latf)
            coordenadas_fin = obtener_coordenadas(punto_fin, loni, lonf, lati, latf)

            resultado = tiempo(coordenadas_inicio[1], coordenadas_inicio[0], coordenadas_fin[1], coordenadas_fin[0])
            if math.isnan(resultado):
                resultado = "nan"
                return (resultado), False
            else:
                
                return (resultado), False
        
        elif comando == "num_puntos":
            n_puntos = num_puntos(lati, loni, latf, lonf)
            return str(n_puntos), False
        
        else:
            return "Comando no reconocido", False
    except Exception as e:
        return f"Error: {str(e)}", False


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Servidor que recibe peticiones de ACO")

    # Definimos argumentos
    parser.add_argument("lati", type=float, help="Latitud inicial del oceano")  
    parser.add_argument("loni", type=float, help="Longitud inicial del oceano")  
    parser.add_argument("latf", type=float, help="Latitud final del oceano")  
    parser.add_argument("lonf", type=float, help="Longitud final del oceano")  
    parser.add_argument("resolucion", type=float, help="Resolucion de los puntos del oceano") 
    parser.add_argument("velocidad", type=float, help="Velocidad del barco")  

    # Parsear los argumentos
    args = parser.parse_args()

    inicio = time.time()

    # Limitar el oceano
    ocean_data.crop_data([args.lati, args.loni, args.latf, args.lonf])

    # Velocidad del barco
    VEL_BARCO = args.velocidad

    # Resolución del mapa
    RESOLUCION = args.resolucion

    # Nombres de los pipes
    pipe_wd11 = "/tmp/pipe_wd11"
    pipe_rd11 = "/tmp/pipe_rd11"

    # Crear los pipes nombrados si no existen
    if not os.path.exists(pipe_wd11):
        os.mkfifo(pipe_wd11)
    if not os.path.exists(pipe_rd11):
        os.mkfifo(pipe_rd11)

    try:
        print("Pipes creados...")
        while True:
            # Leer el mensaje del ACO desde el pipe de entrada
            with open(pipe_wd11, "r") as pipe_in:
                mensaje = pipe_in.read().strip()

                # Validar el formato del mensaje
                if "|" not in mensaje:
                    respuesta = "Error: Formato de mensaje inválido"
                    continuar = False
                else:
                    comando, datos = mensaje.split("|", 1)
                    respuesta, continuar = comunicacion(comando, datos, args.lati, args.loni, args.latf, args.lonf)

                # Escribir la respuesta en el pipe de salida
                with open(pipe_rd11, "w") as pipe_out:
                    pipe_out.write(str(respuesta))

                # Si el comando era "salir", terminar el servidor
                if continuar:
                    print("Cerrando...")
                    fin = time.time()

                    print(f"Tiempo de ejecucion: {fin - inicio} segundos")
                    break

    finally:
        # Eliminar los pipes nombrados al cerrar el servidor
        if os.path.exists(pipe_wd11):
            os.unlink(pipe_wd11)
        if os.path.exists(pipe_rd11):
            os.unlink(pipe_rd11)