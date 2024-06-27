
#!/bin/bash
partype=$1
parname=$2
defaultval=$3
description=$4

infile=$5
parh=$6
parcpp=$7
prev_par=$8

if [[ ! $prev_par ]]; then
	echo "Insufficient input parameters. Values are:
\$1:	partype		C++ data type
\$2:	parname		Name of new parameter
\$3:	defaultval	Default value in C++
\$4:	description	Description for comments
\$5:	infile		path to infile.in
\$6:	parh		path to parameters.h
\$7:	parcpp		path to parameters.cpp
\$8:	prev_par	previous parameter, new lines will be appended after that one"
	exit
fi


## INFILE.IN
echo "INFILE.IN"
if [[ ! $(grep "$parname =" $infile) && $(grep "$prev_par" $infile) ]]; then
	sed -i "/^$prev_par = /a $parname = $defaultval" $infile
	pad=$(printf '%0.1s' " "{1..100})
	padlength=56  ## length of pad string. Currently, the first 56 columns
	auxstr="# $parname = $defaultval"
	auxstr="$auxstr$(printf '%*.*s' 0 $((padlength - ${#parname} - ${#defaultval} - 5 )) "$pad")" ## Need to subtract the length of everything that has already been filled. Including the '# ' and ' = ' in the line (explainging the -5)
	auxstr="$auxstr""## $description "'\['$partype'\]'
	sed -i "/^# $prev_par /a $auxstr" $infile
elif [[ ! $(grep "$prev_par" $infile) ]]; then
	echo " error: the previous parameter $prev_par could not be found in the infile $infile"
else
	echo " error: a parameter by the name $parname is already defined in the infile $infile"
fi

## PARAMETERS.H
echo "PARAMETERS.H"
if [[ ! $(grep "$parname"_ $parh) && $(grep "$prev_par"_ $parh) ]]; then
	lineprev=$(grep -n "$prev_par"_ $parh | head -n1 | cut -f1 -d:)
	assignment_str="$partype $parname""_ = $defaultval;"
	padlength=56
	tabnum=$(( (padlength - ${#assignment_str} - 1) / 4))
	newline=$(echo -e "TAB\t$partype $parname""_ = $defaultval;")
	for i in `seq 1 $tabnum`; do
		newline=$(echo -e "$newline\t")
	done
	newline=$(echo -e "$newline\/\/\/\< $description")
	sed -i "$lineprev"a"$newline" $parh
	sed -i "s/TAB//" $parh

	lineprev=$(grep -n "$prev_par"_ $parh | grep inline | cut -f1 -d:)
	newline=$(echo 'inline '$partype' '$parname'() const { return '$parname'_; };')
	sed -i "$lineprev"a"$newline" $parh
elif [[ ! $(grep "$prev_par"_ $parh) ]]; then
	echo " error: the previous parameter $prev_par could not be found in the parameter header $parh"
else
	echo " error: a parameter by the name $parname is already defined in the parameter header $parh"
fi

## PARAMETERS.CPP
echo "PARAMETERS.CPP"
if [[ ! $(grep "$parname"_ $parcpp) && $(grep "$prev_par"_ $parcpp) ]]; then
	lineprev=$(grep -n "$prev_par""_ = " $parcpp | head -n1 | cut -f1 -d:)
	newline1=$(echo -e 'TAB\t\t} else if (line.rfind("'$parname' =", 0) == 0 ){')
	if [[ "$partype" == "double" || "$partype" == "int" || "$partype" == "bool" ]]; then
		newline2=$(echo -e 'TAB\t\t\t'$parname'_ = std::stod(line.erase(0,line.find("=")+2));')
	elif [[ "$partype" == "std::string" ]]; then
		newline2=$(echo -e 'TAB\t\t\t'$parname'_ = line.erase(0,line.find("=")+2);')
	elif [[ "$partype" == "char" ]]; then
                newline2=$(echo -e 'TAB\t\t\t'$parname'_ = line.at(line.find("=")+2);')
	else
		echo " error: unknown data type $partype. Abort"
		exit
	fi
	sed -i "$lineprev"a"$newline1" $parcpp
	sed -i "$((lineprev+1))"a"$newline2" $parcpp
	sed -i "s/TAB//" $parcpp

	lineprev=$(grep -n "$prev_par"":" $parcpp | grep stdoutfile | head -n1 | cut -f1 -d:)
	newline=$(echo -e 'TAB\tstdoutfile \<\< std::setw(8) \<\< "'$parname':')
	pad=$(printf '%0.1s' " "{1..100})
        padlength=25  # length of pad string
	newline=$(echo -e "$newline$(printf '%*.*s' 0 $((padlength - ${#parname})) "$pad")"'"\t\t\t\<\< '$parname'_')
	padlength=29
	tabnum=$(( ($padlength - ${#parname} + 3 - 1) / 4))  # -1 for the _, +3 for the correct division with rest
        for i in `seq 1 $tabnum`; do
                newline=$(echo -e "$newline\t")
        done    
        newline=$(echo -e "$newline\<\< std::endl;")
	sed -i "$lineprev"a"$newline" $parcpp
        sed -i "s/TAB//" $parcpp

elif [[ ! $(grep "$prev_par"_ $parcpp) ]]; then
        echo " error: the previous parameter $prev_par could not be found in the parameter cpp $parcpp"
else
        echo " error: a parameter by the name $parname is already included in the parameter cpp $parcpp"
fi
