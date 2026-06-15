#```python
import matplotlib.pyplot as plt
import numpy as np
import xml.etree.ElementTree as ET
import os
from scipy.optimize import curve_fit

# --- FUNZIONE DI LETTURA file XML con catalogo indici ---
CATALOGO_INDICI = {
     "N-FK56": 1.43,
     "FK3": 1.46,
     "FK52": 1.48,
     "BK7": 1.51,
     "TIF3": 1.54,
     "LF7": 1.56,
     "BAF4": 1.60,
     "SF12": 1.64,
     "BASF12": 1.67,
     "SF59": 1.95

}

def carica_parametri_automatico():
    path_file = r"C:\PROVE_TESI\DATI_TESI1.optx"
    tree = ET.parse(path_file)
    root = tree.getroot()
    
    #! 1. Legge il primo materiale trovato (riga 22 del file DATI_TESI1.optx, modifica qui il materiale se non si ha quadoa, se lo si ha modificare lì e salvare)
    mat = root.find(".//lens/material")
    nome = mat.get('name') if mat is not None else "Unknown"
    
    #! 2. Legge l'indice dal catalogo associandolo al nome estratto, se si inserisce materiale diverso ad esempio modificando in quadoa usa 1.51
    
    n_mip_0 = CATALOGO_INDICI.get(nome, 1.5)
    #! 3. Legge il passo reticolo direttamente dal file, se non si ha quadoa modifca riga 41 del file DATI_TESI1.optx
    grating = root.find(".//phase[@type='grating']")
    var_lines = grating.find(".//variable[@name='lines_per_mm']")
    lines_per_mm = float(var_lines.get('value')) if var_lines is not None else 1000.0
   

    # Ritorna i 3 valori 
    return nome, lines_per_mm, n_mip_0

# --- ESECUZIONE ---
#restituisce 3 valori per 3 variabili
materiale, d_reticolo, n_mip_0 = carica_parametri_automatico()

wavelengths = np.array([650, 700, 750, 800, 850]) 
#rapporto di variazione dell'indice di rifrazione con la concentrazione, in mL/mg
DN_DC = 0.185 #in mL/mg
# --- PARAMETRI CINETICI ENZIMATICI (GOD) ---
K_M = 21.8  # Costante di Michaelis-Menten per la GOD (espressa in mg/dL per coerenza con l'input)
K_CAT = 100.0  # Costante catalitica (valore indicativo per la simulazione del turnover)
# --- ARRAY STANDARD PER CONFRONTO SCIENTIFICO ---
CONCENTRAZIONI_STANDARD = [50, 60, 70, 80, 90, 140, 150, 160, 170, 180, 190, 200]

#! conversione in mg/dL
# --- INPUT DINAMICO (ARRAY) ---
# --- INPUT DINAMICO ---
q_input = int(input("Quante concentrazioni vuoi simulare? "))
input_str = input(f"Inserisci le {q_input} concentrazioni separate da virgola: ")
lista_conc_utente = [float(c.strip()) for c in input_str.split(',')]
print("Configurazione soglie:")
c_ipo = float(input("Inserisci concentrazione per 'Patologico (es. 60)': "))
c_normale = float(input("Inserisci concentrazione per 'Normoglicemia (es. 90)': "))
c_iper = float(input("Inserisci concentrazione per 'Diabete (es. 140)': "))

# Questa è l'unica lista 'concentrazioni' che devi usare
concentrazioni = sorted(list(set([c_ipo, c_normale, c_iper] + lista_conc_utente)))

input_analita = 10 
# Rimuoviamo C_input_glu e C_effettiva qui, le definiremo solo dove servono




# Creiamo la lista basata solo su quello che decidi tu al momento
#concentrazioni = [c_ipo, c_normale, c_iper]
# --- CALCOLO DINAMICO DELLE SOGLIE ---
# Definisci le tue concentrazioni (es. c_norm, c_iper, c_ipo) prima del grafico
# Assicuriamoci che n_tot sia calcolato per queste soglie
def calcola_angolo_soglia(C_val, n_0, d_ret, dn_dc):
    n_tot = calcola_n_tot(C_val, input_analita, n_0, dn_dc)
    sin_theta = wavelengths / (d_ret * n_tot)
    # Calcoliamo la media degli angoli su tutte le lunghezze d'onda per avere un valore Y unico
    return np.mean(np.degrees(np.arcsin(sin_theta)))
def calcola_n_tot(C_glucosio, C_god, n_base, dn_dc):
    """
    Calcola n_tot applicando la cinetica di Michaelis-Menten.
    C_glucosio e C_god sono in mg/dL.
    """
    V_max = K_CAT * C_god
    # Calcolo del glucosio convertito localmente in gluconolattone
    C_convertita_mg_dL = (V_max * C_glucosio) / (K_M + C_glucosio)
    
    # Conversione in mg/mL per coerenza dimensionale con DN_DC
    C_convertita_mg_mL = C_convertita_mg_dL / 100.0
    
    # Calcolo finale dell'indice di rifrazione
    return n_base + (dn_dc * C_convertita_mg_mL)

# Calcola i valori Y per le linee
y_norm = calcola_angolo_soglia(c_normale, n_mip_0, d_reticolo, DN_DC)
y_iper = calcola_angolo_soglia(c_iper, n_mip_0, d_reticolo, DN_DC)
y_ipo  = calcola_angolo_soglia(c_ipo, n_mip_0, d_reticolo, DN_DC)


C_input_glu_plot = []
angoli_di_diffrazione_plot = []

# Calcolo riferimento fisso (senza concentrazione)
angoli_base = np.degrees(np.arcsin(np.clip(wavelengths / (d_reticolo * n_mip_0), -1, 1)))

for C in concentrazioni:
    n_tot = calcola_n_tot(C, input_analita, n_mip_0, DN_DC)
    sin_theta = wavelengths / (d_reticolo * n_tot)
    angoli = np.degrees(np.arcsin(np.clip(sin_theta, -1, 1)))

    C_input_glu_plot.append(C)
    angoli_di_diffrazione_plot.append(angoli)
    
    pendenza, _ = np.polyfit(wavelengths, angoli, 1)
    
    if C > 0:
        shift_medio = np.mean(angoli_base - angoli)
        sensibilita = shift_medio / C
        print(f"  Shift medio rispetto a base: {shift_medio:.4f}°")
        print(f"  Sensibilità calcolata: {sensibilita:.4f} °/(g/mL)")
    else:
        angoli_base = angoli 
        
    if C in concentrazioni:
        print(f"\n--- ANALISI TARGET: {input_analita} ({C} mg/dL) ---")
        print(f"Angoli diffrazione [°]: {np.round(angoli, 3)}")
        print(f"Pendenza media: {pendenza:.6f} °/nm")
        print(f"Shift medio rispetto a base: {shift_medio:.4f}°")
        print(f"Sensibilità calcolata: {sensibilita:.4f} °/(mg/dL)")

        print(f"  Shift medio rispetto a base: {shift_medio:.4f}°")
        print(f"  Sensibilità calcolata: {sensibilita:.4f} °/(g/mL)")
    else:
            angoli_base = angoli # Memorizziamo per i calcoli successivi
       

print(f"\nAnalisi per: {materiale} | Passo: {d_reticolo:.2f} nm | Indice: {n_mip_0}")

# --- 1. GRAFICO SINGOLO (CURVE DIFFRAZIONE) ---
plt.figure(figsize=(10, 6))
import matplotlib.cm as cm
colors = cm.viridis(np.linspace(0, 0.8, len(lista_conc_utente)))

for idx, C in enumerate(concentrazioni):
    n_tot = calcola_n_tot(C, input_analita, n_mip_0, DN_DC)
    sin_theta = wavelengths / (d_reticolo * n_tot)
    angoli = np.degrees(np.arcsin(np.clip(sin_theta, -1, 1)))
    
    if C in lista_conc_utente:
        # Troviamo l'indice del target per pescare il colore
        target_idx = lista_conc_utente.index(C)
        plt.plot(wavelengths, angoli, label=f'Target: {C} mg/dL', 
                 linewidth=3, color=colors[target_idx], marker='o')
    else:
        # Le linee non-target le lasciamo grigie/alpha per non distrarre
        plt.plot(wavelengths, angoli, label=f'C: {C} mg/dL', alpha=0.3, color='gray')
       

plt.axhline(y=y_norm, color='green', linestyle='--', label=f'Soglia Normale ({c_normale} mg/dL)')
plt.axhline(y=y_iper, color='orange', linestyle='--', label=f'Soglia Iperglicemia ({c_iper} mg/dL)')
plt.axhline(y=y_ipo, color='red', linestyle='--', label=f'Soglia Ipoglicemia ({c_ipo} mg/dL)')
plt.title(f"Sensibilità Angolare: {materiale}")
plt.xlabel("Lunghezza d'onda [nm]")
plt.ylabel("Angolo di diffrazione [°]")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

# --- 2. CALCOLO LINEARITÀ (TERMINALE) ---
lista_c = []
lista_s = []
for C in CONCENTRAZIONI_STANDARD:
    n_tot = calcola_n_tot(C, input_analita, n_mip_0, DN_DC)
    sin_theta = wavelengths / (d_reticolo * n_tot)
    angoli = np.degrees(np.arcsin(sin_theta))
    shift_medio = np.mean(angoli_base - angoli)
    sens = shift_medio / C if C > 0 else 0
    lista_c.append(C)
    lista_s.append(shift_medio)

x = np.array(lista_c)
y = np.array(lista_s)
p = np.polyfit(x, y, 1)
y_pred = np.polyval(p, x)
r_sq = 1 - (np.sum((y - y_pred)**2) / np.sum((y - np.mean(y))**2))

print(f"\n--- VALUTAZIONE TESI: {materiale} ---")
print(f"Coefficiente di linearità (R^2): {r_sq:.4f}")

plt.figure(figsize=(8, 6))

# 1. Disegna i puntini (i dati)
plt.scatter(x, y, color='red', label='Dati Sperimentali', zorder=5)

# 2. Disegna la linea che unisce i puntini (la "piega")
plt.plot(x, y, color='red', linestyle='--', alpha=0.5, label='Andamento sperimentale')

# 3. Disegna la retta di regressione (per far vedere il confronto)
plt.plot(x, y_pred, color='blue', marker='o',  linestyle='-', label=f'Regressione (R^2={r_sq:.4f})')

plt.xlabel('Concentrazione Glucosio [mg/dL]')
plt.ylabel('Shift Angolare [°]')
plt.title(f'Risposta del sensore: {materiale}')
plt.legend()
plt.grid(True)
plt.show()

# --- GRAFICO RIASSUNTIVO  ---
# Dati estratti dai tuoi log
concentrazioni = [65, 90, 180]
shift_nfk56 = [7.18, 7.75, 8.54] # N-FK56
shift_baf4 = [8.08, 8.58, 9.29]  # BAF4

plt.figure(figsize=(8, 6))

# Plot dei dati
plt.plot(concentrazioni, shift_nfk56, marker='o', linestyle='-', linewidth=2, label='N-FK56 ($n=1.43$)')
plt.plot(concentrazioni, shift_baf4, marker='s', linestyle='-', linewidth=2, label='BAF4 ($n=1.60$)')

# Personalizzazione estetica
plt.xlabel('Concentrazione Glucosio [mg/dL]', fontsize=12)
plt.ylabel('Shift Angolare [°]', fontsize=12)
plt.title('Funzione di trasferimento del sensore', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()

# Salva in formato vettoriale (PDF è perfetto per LaTeX)
plt.savefig('figura_y_shift.pdf', bbox_inches='tight')
plt.show()
# --- CONFRONTO MULTI-MATERIALE Predittivo ---
plt.figure(figsize=(10, 6))

# Definiamo le concentrazioni per cui vogliamo il confronto
concentrazioni_confronto = [50, 70, 90, 100, 126, 140, 200]

for nome_mat, n_0 in CATALOGO_INDICI.items():
    sens_list = []
    for C in concentrazioni_confronto:
        # Ricalcoliamo la sensibilità per questo specifico materiale
        n_tot = calcola_n_tot(C, input_analita, n_mip_0, DN_DC)
        sin_theta = wavelengths / (d_reticolo * n_tot)
        angoli = np.degrees(np.arcsin(sin_theta))
        angoli_base = np.degrees(np.arcsin(wavelengths / (d_reticolo * n_0)))
        
        shift_medio = np.mean(angoli_base - angoli)
        sens = shift_medio / C if C > 0 else 0
        sens_list.append(sens)
    
    # Plottiamo la curva per questo materiale
    plt.plot(concentrazioni_confronto, sens_list, marker='o', label=nome_mat)

plt.title("Confronto Sensibilità Materiali")
plt.xlabel("Concentrazione [mg/dL]")
plt.ylabel("Sensibilità [°/(mg/dL)]")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
# --- SALVATAGGIO CORRETTO ---
percorso_log = r"C:\PROVE_TESI\log_risultati.txt"
with open(percorso_log, "a") as f:
    f.write(f"\n--- RUN: {materiale} | Passo: {d_reticolo:.2f} nm ---\n")
    f.write(f"Analita: {input_analita} | DN/DC: {DN_DC} mL/mg\n")

    for C in concentrazioni:
        # Calcoli ricalcolati correttamente per il log
        n_tot = calcola_n_tot(C, input_analita, n_mip_0, DN_DC)
        sin_theta = wavelengths / (d_reticolo * n_tot)
        angoli = np.degrees(np.arcsin(sin_theta))
        pendenza, _ = np.polyfit(wavelengths, angoli, 1)
        shift_medio = np.mean(angoli_base - angoli)
        sensibilita = shift_medio / C if C > 0 else 0
        r_sq = 1 - (np.sum((y - y_pred)**2) / np.sum((y - np.mean(y))**2))
        nome = materiale #usato per poter scrivere l'indice nel log
        n_mip_0 = CATALOGO_INDICI.get(nome, 1.5)
        # SCRITTURA PULITA (Niente input qui dentro!)
        f.write(f"indice di rifrazione: {n_mip_0}| Conc: {C} mg/dL|Array di concentrazioni: {input_str} Pendenza: {pendenza:.6f} | Angoli: {np.round(angoli, 3)} | Shift: {shift_medio:.4f} | Sensibilità: {sensibilita:.6f}\n| Linearità (R^2): {r_sq:.4f}\n")
    
    f.write("------------------------------------------\n")
    f.flush() # Forza lo svuotamento del buffer
    os.fsync(f.fileno()) # Forza la scrittura fisica sul disco



materiali = ['N-FK56', 'FK3', 'FK52', 'BK7', 'BAF4']

# Sensibilità estratte dai log per la concentrazione target di ogni run
sens_ipo = [0.1105, 0.1131, 0.1148, 0.1172, 0.1244]    # Target 58
sens_norm = [0.0790, 0.0805, 0.0815, 0.0830, 0.0872]   # Target 100
sens_diab = [0.0612, 0.0623, 0.0630, 0.0640, 0.0670]   # Target 135

plt.figure(figsize=(10, 6))
plt.plot(materiali, sens_ipo, 'o-', label='Ipoglicemia (58 mg/dL)', color='red')
plt.plot(materiali, sens_norm, 's-', label='Normoglicemia (100 mg/dL)', color='green')
plt.plot(materiali, sens_diab, 'd-', label='Diabete (135 mg/dL)', color='blue')

plt.title('Sensibilità del sistema al variare del substrato e della condizione clinica')
plt.xlabel('Materiale Substrato')
plt.ylabel('Sensibilità [°/(mg/dL)]')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()



print("Salvataggio completato! Controlla il file .txt sul Desktop.")
#```