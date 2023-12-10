for ((C = 1; C < 10001; C *=10)); do
    for ((K=2 ; K < 11; K+=1)); do
        for ((HMAX=2 ; HMAX < 11; HMAX+=2)); do
            hmax=$HMAX
            primero="Archivos/k"
            k=$K
            segundo="_c"
            c=$(awk "BEGIN { printf \"%.2f\", $C / 100 }")
            underscore="_"
            tercero="C.cpp"
            cuarto="C.dat"
            quinto=""
            sexto=""

            nombre=$primero$k$segundo$c$underscore$hmax$tercero
            nombreData="k"$k$segundo$C$underscore$hmax$cuarto
            nombreK=$primero$quinto$k$segundo$C$underscore$hmax$tercero
            nombreC=$primero$k$segundo$sexto$c$underscore$hmax$tercero

            echo $nombre

            # Create temporary files for intermediate results
            temp_file1=$(mktemp)
            temp_file2=$(mktemp)
            temp_file3=$(mktemp)
            temp_file4=$(mktemp)

            sed "48s/.*/#define HMAX    $hmax /" k2_c001_2C.cpp > "$temp_file1"
            sed "49s/.*/#define K       $k /" "$temp_file1" > "$temp_file2"
            sed "52s/.*/#define c	    $c /" "$temp_file2" > "$temp_file3"
            sed "56s/.*/#define namearch        $nombreData /" "$temp_file3" > "$temp_file4"
            # Perform other sed operations similarly, using temporary files

            # Rename the final output file
            mv "$temp_file4" "$nombreK"

            # Clean up temporary files
            rm "$temp_file1"
            rm "$temp_file2"
            rm "$temp_file3"
        done    
    done
done