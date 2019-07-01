#!/bin/bash
# create scheme file 

Usage_exit () {

echo "Usage : fsl2scheme.sh <bvals> <bvecs> <bscale> <output> [-r]"
echo "Option -r : round bval at the nearest hundredth"
exit 1

}

formbvec () {

if [ "`echo $1 | head -c 1`" = "-" ] ; then
  echo $1 | awk '{printf "  %.6f",$1}'
else
  echo $1 | awk '{printf "   %.6f",$1}'
fi

}

formbval () {

if [ $round != 1 ] ; then
 echo $1 | awk '{printf "   %.3e",$1*'$bscale'}' | sed -e 's/e+/E/'
else
 i=`echo $1 | awk '{printf "%d",$1}'`
 echo "(($i + 50)/100) * 100" | bc | awk '{printf "   %.3e",$1*'$bscale'}' | sed -e 's/e+/E/'
fi

}

[ "$3" = "" ] && Usage_exit

bval=(`cat $1 | awk 'NR==1 {print}'`)
bvecx=(`cat $2 | awk 'NR==1 {print}'`)
bvecy=(`cat $2 | awk 'NR==2 {print}'`)
bvecz=(`cat $2 | awk 'NR==3 {print}'`)
bscale=$3
output=$4

if [ "$5" = "-r" ] ; then
  round=1
fi

num=`echo ${bval[@]} | wc -w`

if [ -e $output ] ; then rm $output ;fi
echo -e "# Scheme file created by fsl2scheme.sh\nVERSION: BVECTOR" > $output
i=0
while [ $i -lt $num ] ; do
  bvecxf=`formbvec ${bvecx[$i]}`
  bvecyf=`formbvec ${bvecy[$i]}`
  bveczf=`formbvec ${bvecz[$i]}`
  bvalf=`formbval ${bval[$i]}`
  echo "${bvecxf}${bvecyf}${bveczf}${bvalf}" >> $output
  i=`expr $i + 1`
done
