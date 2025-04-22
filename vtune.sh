MANNA_BINARY="./tiny_manna"

start=1024              # 1K
end=$(( 1024*1024 ))    # 1M
n=$start

while [ $n -le $end ]
do
	size=$(($n*4/1024))
        echo "Tamaño " $size " KBytes. N vale: " $n
        make clean
        make tiny_manna N=$n CXX=clang++

        #execute
        vtune -collect hpc-performance -r report$size $MANNA_BINARY

        n=$(($n*2))
done