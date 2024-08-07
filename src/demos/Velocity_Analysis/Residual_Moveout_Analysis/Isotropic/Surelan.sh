#! /bin/sh
# Residual moveout analyses for the cip gathers
# residual moveout is represented by z(h)^2 = z(0)^2+r*h^2
# 	where h is half an offset
# Authors: Zhenyue
# NOTE: Comment lines preceeding user input start with  ##
#set -x

## Set parameters
input=kd.data.su
rpicks=res.p1
cdpmin=800
cdpmax=2000
dcdp=600
fold=5

## Set r-parameter sampling
nr=51 dr=0.01 fr=-0.25


## Do the residual moveout analyses.
echo 
echo "Pick r-parameters by moving mouse and typing 's', type 'Q' when done"
echo 
echo "Use the plot from Xmig as a guide to picking"

cdp=$cdpmin
while [ $cdp -le $cdpmax ]
do
	ok=false
	while [ $ok = false ]
	do
		echo "Starting Residual moveout analysis for cip $cdp"
		suwind <$input key=cdp min=$cdp max=$cdp count=$fold |
		suxwigb title="cdp=$cdp" xbox=10 label1="time (s) " \
			label2="offset" key=offset &
		suwind <$input key=cdp min=$cdp max=$cdp count=$fold |
		surelan nr=$nr dr=$dr fr=$fr |
		suximage bclip=0.95 wclip=0.0 f2=$fr d2=$dr cmap=hsv2 \
			label1="Depth (m)" label2="r-parameter " xbox=600 \
			title="r-parameter Scan for CIP $cdp" \
			mpicks=mpicks.$cdp blank=.7

		pause

		/bin/echo -n "Picks OK? (y/n) " >/dev/tty
		read response
		case $response in
		n*) ok=false ;;
		*) ok=true ;;
		esac

	done </dev/tty
	cdp=`bc -l <<END
		$cdp + $dcdp
END`

done

set +x


### Combine the individual picks into a composite par file
echo "Editing pick files ..."
>$rpicks
cdp=$cdpmin
while [ $cdp -le $cdpmax ]
do
	/bin/echo -n "cip=$cdp," >temp
	cat temp mpicks.$cdp>>$rpicks
	cdp=`bc -l <<END
		$cdp + $dcdp
END`
done

sed "s/  /,/g"<$rpicks >cig.par

echo "dzdv par file: cig.par is ready"


### Clean up
cdp=$cdpmin
while [ $cdp -le $cdpmax ]
do
	rm mpicks.$cdp 
	cdp=`bc -l <<END
		$cdp + $dcdp
END`
done
rm temp $rpicks
