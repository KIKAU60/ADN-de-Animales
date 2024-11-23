import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from Bio.Seq import Seq

# Función para generar la doble hélice
def generar_helice_adn(secuencia_adn):
    """
    Genera una visualización 3D de una doble hélice de ADN a partir de una secuencia.
    """
    colores = {'A': 'blue', 'T': 'red', 'C': 'green', 'G': 'yellow'}

    # Número de pasos o puntos por vuelta de la hélice
    pasos_por_vuelta = 10
    vueltas = len(secuencia_adn) // pasos_por_vuelta

    # Generar coordenadas en 3D para la doble hélice
    t = np.linspace(0, 4 * np.pi, len(secuencia_adn))
    x = np.sin(t)
    y = np.cos(t)
    z = np.linspace(0, 1, len(secuencia_adn))

    # Crear la visualización
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    for i, base in enumerate(secuencia_adn):
        ax.scatter(x[i], y[i], z[i], color=colores[base], s=100, label=base if i == 0 else "")

    # Dibujar las conexiones entre las dos cadenas
    for i in range(0, len(secuencia_adn) - 1, 2):
        ax.plot([x[i], x[i+1]], [y[i], y[i+1]], [z[i], z[i+1]], color='black', lw=1)

    # Personalización
    ax.set_title("Representación 3D de la Doble Hélice de ADN")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.grid(False)

    # Ajustes de la vista
    ax.view_init(30, 60)
    plt.legend()

    # Mostrar la visualización
    plt.show()

# Función para obtener los codones de una secuencia de ADN
def obtener_codones(secuencia_adn):
    """
    Extrae los codones de una secuencia de ADN.
    """
    return [secuencia_adn[i:i+3] for i in range(0, len(secuencia_adn), 3)]

# Función para graficar los codones
def graficar_codones(codones):
    """
    Muestra la frecuencia de los codones en un gráfico de barras.
    """
    frecuencia_codones = Counter(codones)
    codones, frecuencias = zip(*sorted(frecuencia_codones.items()))

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(codones, frecuencias, color='purple')

    ax.set_title('Frecuencia de Codones')
    ax.set_xlabel('Codón')
    ax.set_ylabel('Frecuencia')
    plt.xticks(rotation=90)
    plt.show()

# Función para calcular la proporción de nucleótidos
def calcular_proporcion_nucleotidos(secuencia_adn):
    """
    Calcula y muestra la proporción de nucleótidos en la secuencia de ADN.
    """
    secuencia = Seq(secuencia_adn)

    count_a = secuencia.count('A')
    count_t = secuencia.count('T')
    count_c = secuencia.count('C')
    count_g = secuencia.count('G')

    total_nucleotidos = len(secuencia)

    proporcion_a = count_a / total_nucleotidos
    proporcion_t = count_t / total_nucleotidos
    proporcion_c = count_c / total_nucleotidos
    proporcion_g = count_g / total_nucleotidos

    nucleotidos = ['Adenina (A)', 'Timina (T)', 'Citosina (C)', 'Guanina (G)']
    proporciones = [proporcion_a, proporcion_t, proporcion_c, proporcion_g]
    colores = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    plt.figure(figsize=(8, 8))
    plt.pie(proporciones, labels=nucleotidos, autopct='%1.1f%%', startangle=140, colors=colores)
    plt.title("Proporción de Nucleótidos en la Secuencia de ADN")
    plt.show()

# Función principal para ejecutar todas las visualizaciones
def main():
    # Solicitar la secuencia de ADN
    secuencia_adn = input("Introduce la secuencia de ADN del animal: ")

    # Validar secuencia (solo A, T, C, G)
    if not all(base in 'ATCG' for base in secuencia_adn):
        print("La secuencia de ADN contiene caracteres no válidos. Asegúrate de que solo contenga A, T, C y G.")
        return

    # Generar y visualizar la doble hélice
    generar_helice_adn(secuencia_adn)

    # Obtener y graficar los codones
    codones = obtener_codones(secuencia_adn)
    graficar_codones(codones)

    # Calcular y mostrar la proporción de nucleótidos
    calcular_proporcion_nucleotidos(secuencia_adn)

# Ejecutar el programa
if __name__ == "__main__":
    main()
