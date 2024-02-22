# README

## Installazione ambiente
Aprire con un editor il file `setup_env_satellite.sh` e sostituire il path assoluto corretto della working directory (`fp_cwd`). 
Aprire un terminale nella cartella in cui si trova il file `setup_env_satellite.sh`. Dare permessi di esecuzione con la linea di comando
```
chmod +x ./setup_env_satellite.sh
```
ed eseguire il file con
```
./setup_env_satellite.sh
```
L'eseguibile installa `miniconda`, che è un package manager per Python, e crea un ambiente virtuale e locale (all'interno della cartella `conda` nella working directory).

Per attivare l'ambiente virtuale e installare tutti i pacchetti necessari a far girare i codici spostarsi nella cartella `conda` e lanciare il file di settings:
```
cd conda
source satellite_settings
```
Basta usare il file di settings solo yuna volta, poi successivamente si può attivare l'ambiente con il comando:
```
conda activate satellite_libraries
```
A questo punto in linea di comando si evidenzia (satellite_libraries) che è l'ambiente in cui si vuole lavorare.
Per disattivarlo:
```
conda deactivate
```

Nota: per installare pacchetti o far girare i codici, è necessario attivare l'ambiente da terminale ed eseguire le operazioni su quel terminale.
Attenzione: meglio non usare il file yml per creare l'ambiente: ci vuole troppo tempo. I pacchetti sono hardcoded nel file di setup. 


***

## Run codici
L'ambiente contiene già tutti i pacchetti necessari a far girare i codici.
Attivare l'ambiente.
Attivare Jupyter Lab nella working directory,
```
cd */codes_repietra
cd conda
source satellite_settings
cd ..
jupyter-lab
```
A questo punto si apre una pagina del browser con tutta la working directory e l'editor dei notebook (i codici interattivi che possono girare su Jupyter Lab). Basta andare nella cartella "Codici" e aprire il notebook `main.ipynb`. Prima di farlo girare, assicurarsi di sistemare i path giusti delle cartelle nel file di configurazione `configuration.json` e i flag di salvataggio (vedi sezione configuration.json).

Per far girare il codice, nel menu in alto Run -> Run all cells.
Ne caso in cui si facciano delle modifiche al file di configurazione o altri file importati dal codice, bisogna resettare il kernel per ricaricare i moduli: Kernel -> Restart kernel and run all cells. 


## configuration.json
Il file di configurazione deve essere impostato con tutte le informazioni sui path e sui parametri che si vuole calibrare.
Per aprirlo correttamente, usare un editor di testo in locale (da pc) oppure su Jupyter Lab right-click -> Open with -> Editor (JSON viene fuori solo in modalità visualizzazione).
Tutte le entrate del file di configurazione vengono lette dal notebook main.ipynb in una delle prime celle. Se si vuole modificare il file di configurazione (e.g. elimina, aggiungi una variabile, etc) bisogna accertarsi di cambiare il corrispondente nel notebook. Allo stato originale il file di configurazione risponde alla maggior parte delle esigenze che si possono avere nell'utilizzo del codice.

Variabili da modificare:
- tutte quelle nei campi `options` e i path in `paths`;
- `opt_save` ("True"/"False") per gli output, conferma il salvataggio di plot e tabelle;
- `automate` ("True"/"False") alla fine del file, fa girare la procedura e i salvataggi in modo automatico.

Nota: non è necessario settare un filename per i plot o per le tabelle, in quanto il setting è automatizzato e impostato dall'impostazione `filename_template`, che può essere modificata aggiungendo tra graffe i nomi di altre variabili presenti nel file di configurazione.
Alla fine di ogni file name viene apposto "add_description".
