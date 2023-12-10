#sed 'Ns/.*/replacement-line/' file.txt > new_file.txt
#Líneas: 48, 49, 52, 56
#define HMAX    2    	// numero maximo de capas que maneja el programa
#define K       2      	// tama�o del kmero
#define c	    0.01   					// cociente q1/qi
#define namearch        "k2_c001_2C.dat"
#0.01 0.1 1 10 1000
#1/100 10/100 100/100 1000/100 10000/100

for ((HMAX=2 ; HMAX < 11; HMAX+=2)); do
    for ((K=2 ; K < 11; K+=1)); do
        for ((C = 1; C < 10001; C *=10)); do
            hmax=$HMAX
            primero="Archivos/k"
            k=$K
            segundo="_c"
            c=$C
            underscore="_"
            tercero="C.cpp"
            cuarto="C.dat"

            nombre=$primero$k$segundo$c$underscore$hmax$tercero
            nombreData=$primero$k$segundo$c$underscore$hmax$cuarto
            echo $nombre
            sed "48s/.*/#define HMAX    $hmax /" k2_c001_2C.cpp > $nombre
            sed -i "49s/.*/#define K       $k /" $nombre
            #sed "52s/.*/#define c	    $c /" k2_c001_2C.cpp > $nombre
            #sed "56s/.*/#define namearch        "'$nombreData'" /" k2_c001_2C.cpp> $nombre
        done    
    done
done
