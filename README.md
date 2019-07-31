# IMAT
Iterative Methods Analysis Tool

Tool Matlab per la risoluzione di sistemi lineari con i metodi iterativi di Jacobi, Gauss-Seidel e il gradiente (coniugato e/o precondizionato)

Lo strumento di analisi è stato progettato allo scopo di essere di semplice e rapido utilizzo ed inoltre particolarmente flessibile e personalizzabile. All’avvio del programma, l’utente potrà:

1. selezionare la tipologia di matrice che desidera utilizzare per i test

2. definire il numero di test (step) di diversa dimensione che si vuole considerare (es. 10, 100, 1000)

3. modificare i parametri di generazione della matrice predefiniti (es. condizionamento e sparsità)

4. configurare i parametri per l’esecuzione del test (es. tolleranza e numero massimo di iterazioni)

5. visualizzare in formato tabellare i risultati relativi a tempi, iterazioni ed errori relativi

6. visualizzare infine le seguenti **tipologie di grafici**:
   * grafici in scala logaritmica per i residui
   * grafici a barre per i tempi e le iterazioni
   * grafici a barre per gli errori relativi
 
 ## Funzioni
 Le funzioni implementate sono le seguenti:
  * **IterMethodsTool**: funzione semi-automatizzata per:
      * acquisire l’input utente per la configurazione del problema
      * effettuare le chiamate alle singole funzioni dei metodi iterativi
      * memorizzare e visualizzare i risultati
  * **MatrixCreator**: funzione per la creazione delle matrici dalla differente struttura
  * **Varie**: funzioni specifiche per i singoli metodi iterativi

## Installazione
Per utilizzare il tool è sufficiente eseguire il file **IterMethodTool.m** con Matlab. Dopo aver inserito i parametri necessari al calcolo dei diversi metodi, i risultati verranno visualizzati su schermo con dei grafici.

## Test
Il tool è stato provato con Matlab _R2019a_ senza riscontrare alcun problema.

------------------------------------

# IMAT (english)
Iterative Methods Analysis Tool

Tool Matlab for the resolution of linear systems with the iterative methods of Jacobi, Gauss-Seidel and the gradient (conjugate, and/or preconditioned)

The analysis tool has been designed in order to be simple and quick to use and also particularly flexible and customizable. When the program starts, the user can:

1. select the type of matrix to use for the tests

2. define the number of tests (steps) of different dimensions to be considered (eg 10, 100, 1000)

3. modify the predefined matrix generation parameters (eg conditioning and sparsity)

4. configure the parameters for the execution of the test (eg tolerance and maximum number of iterations)

5. display the results relating to times, iterations and relative errors in tabular format

6. finally display the following **types of charts**:
    * graphs in logarithmic scale for residuals
    * bar graphs for times and iterations
    * bar graphs for relative errors
 
 ## Functions
 The functions implemented are the following:
 * **IterMethodsTool**: semi-automated function for:
      * acquire the user input for the configuration of the problem
      * make calls to the individual functions of the iterative methods
      * memorize and display the results
 * **MatrixCreator**: function for creating matrices from different structure
 * **Miscellaneous**: specific functions for individual iterative methods

## Installation
To use the tool, simply run the **IterMethodTool.m** file with Matlab. After entering the parameters necessary for the calculation of the different methods, the results will be displayed on the screen with graphs.

## Test
The tool was tested with Matlab _R2019a_ without encountering any problems.
