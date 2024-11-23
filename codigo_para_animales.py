import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from Bio.Seq import Seq

# Importar el archivo con las secuencias de ADN
from secuencias_adn import secuencias_adn

# Función para generar la doble hélice con la estructura secundaria
def generar_helice_adn(secuencia_adn):
    """
    Genera una visualización 3D de la doble hélice de ADN con la estructura secundaria.
    """
    colores = {'A': 'blue', 'T': 'red', 'C': 'green', 'G': 'yellow'}
    
    # Número de pasos o puntos por vuelta de la hélice
    pasos_por_vuelta = 10
    vueltas = len(secuencia_adn) // pasos_por_vuelta
    
    # Generar coordenadas en 3D para la doble hélice
    t = np.linspace(0, 4 * np.pi, len(secuencia_adn))  # Variable para la espiral
    x1 = np.sin(t)  # Coordenadas para la primera cadena
    y1 = np.cos(t)
    z1 = np.linspace(0, 1, len(secuencia_adn))

    x2 = np.sin(t + np.pi)  # Coordenadas para la segunda cadena (desfasada 180 grados)
    y2 = np.cos(t + np.pi)
