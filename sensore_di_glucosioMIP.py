#```python
import matplotlib.pyplot as plt
import numpy as np
import xml.etree.ElementTree as ET
import os

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


#! conversione in mg/dL
C_input_glu = float(input("Inserisci la concentrazione di zucchero da simulare in mg/dL  (es. 70, 126, 140):")) # valori di glicemia a digiuno da testare, inserisci 70 per ièpoglicemia, 126 per prediabete, 140 per diabete
C_effettiva = C_input_glu /100 # conviene comunque lavorare in mg/mL, ma l'interfaccia sarà in mg/dL per comodità dell'utente
input_analita = 10 #! Fissato 

# --- INPUT DINAMICO ---
print("Configurazione simulazione:") #? queste sono le soglie cliniche, verranno usate èper delimitare il comportamento del sensore entro questi range
c_ipo = float(input("Inserisci concentrazione per 'Patologico' (es. 60): "))
c_normale = float(input("Inserisci concentrazione per 'Normoglicemia' (es. 100): "))
c_iper = float(input("Inserisci concentrazione per 'Iperglicemia' (es. 200): "))


# Creiamo la lista basata solo su quello che decidi tu al momento
concentrazioni = [c_ipo, c_normale, c_iper]
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




if C_input_glu not in concentrazioni:
    concentrazioni.append(C_input_glu)
#liste per salvare le concentrazioni in mg/dL e i corrispondenti angoli di diffrazione, per poi salvarli 
C_input_glu_plot = []
angoli_di_diffrazione_plot = []


angoli_base = np.degrees(np.arcsin(wavelengths / (d_reticolo * n_mip_0))) # Angoli di base senza concentrazione, per calcolare lo shift

for C in concentrazioni:
    # 1. Converti in unità interna (mg/mL)
    c_interna = C_effettiva
    n_tot = calcola_n_tot(C, input_analita, n_mip_0, DN_DC) # facciamo in modo che C sia in mg/mL e DN_DC in mL/mg, quindi moltiplichiamo C per 1/100 per ottenere la concentrazione in mg/dL
    sin_theta = wavelengths / (d_reticolo * n_tot)
    angoli = np.degrees(np.arcsin(sin_theta))
    angoli_base = np.degrees(np.arcsin(wavelengths / (d_reticolo * n_mip_0))) # Angoli di base senza concentrazione, per calcolare lo shift

    # Salvataggio nelle liste per il grafico
    C_input_glu_plot.append(C_input_glu)
    angoli_di_diffrazione_plot.append(angoli)
    # Calcolo pendenza (Dispersion Rate)
    pendenza, _ = np.polyfit(wavelengths, angoli, 1)
    # Calcolo Sensibilità rispetto a C=0
    if C > 0:
            shift_medio = np.mean(angoli_base - angoli)
            sensibilita = shift_medio / C
            print(f"  Shift medio rispetto a base: {shift_medio:.4f}°")
            print(f"  Sensibilità calcolata: {sensibilita:.4f} °/(g/mL)")
    else:
            angoli_base = angoli # Memorizziamo per i calcoli successivi
        
        # Stampa dettagliata nel terminale
    if abs(C - C_input_glu) < 0.01:
            print(f"\n--- ANALISI TARGET: {input_analita} ({C} mg/dL) ---")
            print(f"Angoli diffrazione [°]: {np.round(angoli, 3)}")
            print(f"Pendenza media: {pendenza:.6f} °/nm")
            print(f"Shift medio rispetto a base: {shift_medio:.4f}°")
            print(f"Sensibilità calcolata: {sensibilita:.4f} °/(mg/dL)")
    

print(f"\nAnalisi per: {materiale} | Passo: {d_reticolo:.2f} nm | Indice: {n_mip_0}")

# --- 1. GRAFICO SINGOLO (CURVE DIFFRAZIONE) ---
plt.figure(figsize=(10, 6))
for C in concentrazioni:
    n_tot = calcola_n_tot(C, input_analita, n_mip_0, DN_DC)
    sin_theta = wavelengths / (d_reticolo * n_tot)
    angoli = np.degrees(np.arcsin(sin_theta))
    if C == C_input_glu:
        plt.plot(wavelengths, angoli, label=f'Target: {C} mg/dL', linewidth=3, color='black', marker='o')
    else:
        plt.plot(wavelengths, angoli, label=f'C: {C} mg/dL', alpha=0.6)

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
for C in concentrazioni:
    n_tot = calcola_n_tot(C, input_analita, n_mip_0, DN_DC)
    sin_theta = wavelengths / (d_reticolo * n_tot)
    angoli = np.degrees(np.arcsin(sin_theta))
    shift_medio = np.mean(angoli_base - angoli)
    sens = shift_medio / C if C > 0 else 0
    lista_c.append(C)
    lista_s.append(sens)

x = np.array(lista_c)
y = np.array(lista_s)
p = np.polyfit(x, y, 1)
y_pred = np.polyval(p, x)
r_sq = 1 - (np.sum((y - y_pred)**2) / np.sum((y - np.mean(y))**2))

print(f"\n--- VALUTAZIONE TESI: {materiale} ---")
print(f"Coefficiente di linearità (R^2): {r_sq:.4f}")

# --- GRAFICO RIASSUNTIVO  ---
plt.figure(figsize=(10, 6))
# Usiamo i dati che hai già calcolato nella sezione precedente
plt.scatter(lista_c, lista_s, color='blue', marker='o', label='Punti di test')
plt.plot(lista_c, y_pred, color='red', linestyle='--', label=f'Regressione (R^2={r_sq:.4f})')

plt.title(f"Linearità Sensibilità: {materiale}")
plt.xlabel("Concentrazione [mg/dL]")
plt.ylabel("Sensibilità [°/(mg/dL)]")
plt.legend()
plt.grid(True)
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
        # SCRITTURA PULITA (Niente input qui dentro!)
        f.write(f"Conc: {C} mg/dL| Concentrazione di glucosio: {C_input_glu} | Pendenza: {pendenza:.6f} | Angoli: {np.round(angoli, 3)} | Shift: {shift_medio:.4f} | Sensibilità: {sensibilita:.6f}\n")
    
    f.write("------------------------------------------\n")
    f.flush() # Forza lo svuotamento del buffer
    os.fsync(f.fileno()) # Forza la scrittura fisica sul disco

print("Salvataggio completato! Controlla il file .txt sul Desktop.")
#```