# README


**Analisi dell’espressione di circRNA umani e virali nei precursori neuronali infettati da virus Zika e Citomegalovirus umano**

- **circRNA_analysis**
  - **illumina_analysis**
    - *data_preprocessing*
      - run_illumina_fastp.sh
    - *alignment*
      - CIRI2+CIRI-AS+CIRI-full.sh
    - *read_counting*
      - count_reads_illumina.py
      - tot_circRNA_illumina.py
    - *output/*
  - **nanopore_analysis**
    - *data_preprocessing*
    - *alignment*
    - *circRNA_detection*
      - long_read_circRNA_run.sh
      - formatted.py
    - *read_counting*
      - fasta_intermedi.py
    - *normalization*
      - normalize_expression.py
    - *output/*
  - README.md



**Autore**: Martina Zambon  
**Relatore**: Prof. Enrico Lavezzo  
**Anno Accademico**: 2023-2024

---

## Descrizione

Questo progetto di tesi ha lo scopo di analizzare i profili di espressione dei circRNA in cellule precursori neuronali infettate da Zika Virus (ZIKV) e Citomegalovirus umano (CMV) a diverse condizioni temporali (24h, 48h, 72h). Lo studio utilizza tecniche di sequenziamento sia a lettura breve (Illumina) che a lettura lunga (Nanopore) per confrontare e integrare i risultati ottenuti da entrambe le piattaforme. 

I dati ottenuti sono stati analizzati utilizzando una pipeline bioinformatica appositamente costruita, che ha permesso di identificare e analizzare l’espressione differenziale dei circRNA. I risultati forniscono insight sui meccanismi molecolari coinvolti nella risposta delle cellule all'infezione virale e potrebbero avere implicazioni importanti nello sviluppo di biomarcatori per malattie congenite.

---

## Struttura della Pipeline

### 1. **Pre-processing dei Dati di Sequenziamento**

Le letture ottenute dalle piattaforme di sequenziamento sono state preprocessate per garantire una qualità adeguata prima dell'analisi. Questo include:
- **Filtraggio e trimming** delle letture grezze usando [fastp](https://github.com/OpenGene/fastp) per la rimozione di basi di bassa qualità e adattatori.
- **Allineamento** delle letture al genoma di riferimento umano (GRCh38) usando [BWA-MEM](http://bio-bwa.sourceforge.net/) per i dati Illumina, e [pblat](https://icebert.github.io/pblat/) per i dati Nanopore.
- **Rimozione di duplicati** e pulizia dei dati per migliorare la precisione delle analisi downstream.

### 2. **Identificazione dei circRNA**
Per l’identificazione dei circRNA sono stati usati i seguenti strumenti:
- **[CIRI2](https://github.com/bioinfo-biols/CIRI-cookbook/blob/master/source/CIRI2.md)** per i dati Illumina, che identifica giunzioni di back-splicing dai dati di lettura breve.
- **[long_read_circRNA](https://github.com/omiics-dk/long_read_circRNA)** per i dati Nanopore, specificamente progettato per catturare circRNA a lettura lunga.

### 3. **Analisi dell’Espressione Differenziale**
L’analisi dell’espressione differenziale è stata eseguita utilizzando un approccio statistico rigoroso per determinare i circRNA significativamente espressi in modo differenziale tra cellule infettate da CMV/ZIKV e controlli:
- Normalizzazione dei dati (TPM) per tener conto delle differenze nella profondità di sequenziamento.
- Uso di **modelli di distribuzione di Poisson negativa** per l’analisi statistica dell’espressione differenziale.
- Visualizzazione dei risultati con **volcano plots** e **heatmaps** per evidenziare i circRNA con variazioni significative nell’espressione.

### 4. **Visualizzazione e Interpretazione dei Risultati**
Sono stati generati grafici come:
- **Volcano plots** per mostrare l'espressione differenziale dei circRNA.
- **Heatmaps** per visualizzare i pattern di espressione.
- **Boxplots** per confrontare i campioni nelle diverse condizioni temporali.
  
---

## Requisiti Tecnici

Per eseguire la pipeline, sono necessari i seguenti strumenti e pacchetti:

- **Python 3.x**
  - pandas
  - matplotlib
  - seaborn
- **R** con i seguenti pacchetti:
  - DESeq2
  - ggplot2
- **Bash** (per gli script di preprocessing e allineamento)
- **fastp** per il trimming e la pulizia dei dati
- **BWA-MEM** per l’allineamento delle letture Illumina
- **pblat** per l’allineamento delle letture Nanopore
- **CIRI2** per l’identificazione di circRNA in letture Illumina
- **long_read_circRNA** per l’identificazione di circRNA in letture Nanopore

---

## Dati di Input e Output

- **Input**: I file di input sono costituiti dalle letture di sequenziamento preprocessate (.fastq) e dagli allineamenti (.bam).
- **Output**: Gli output della pipeline includono:
  - File di identificazione dei circRNA (`.txt`).
  - Volcano plots e heatmaps delle analisi di espressione differenziale.
  - Tabelle dei circRNA significativamente espressi in modo differenziale.

---

## Contatti

 
**Martina Zambon**  
Email: martina.zambon.7@studenti.unipd.it
