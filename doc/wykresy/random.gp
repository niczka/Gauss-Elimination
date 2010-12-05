set terminal postscript enhanced 
set output 'random_plot.eps'
set key inside right bottom vertical Right noreverse enhanced 
set samples 10000
set xlabel "rozmiar macierzy"
set ylabel "||b-Ax||"
set grid
set title "Eliminacja Gaussa - macierz losowa"
set offsets 0
set multiplot
plot [0:100][0:0.000001] 'random_without' ti "Bez wyboru" with lines,\
'random_from_row' ti "Wybor z wiersza" with lines,\
'random_from_col' ti "Wybor z kolumny" with lines,\
'random_full' ti "Pelny wybor" with lines
