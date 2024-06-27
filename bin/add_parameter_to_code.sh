
#!/bin/bash
type=$1
parname=$2
defaultval=$3
description=$4
if [[ ! $description ]]; then
	echo "Insufficient input parameters. Values are:
\$1:	type		C++ data type
\$2:	parname		Name of new parameter
\$3:	defaultval	Default value in C++
\$4:	description	Description for comments"
	exit
fi
echo "IN DEFAULTS"
echo "$parname = $defaultval"
echo "# $parname = $defaultval				## $description"
echo ""

echo "IN H"
echo "$type $parname"_" = $defaultval; ///< $description"
echo 'inline '$type' '$parname'() const { return '$parname'_; };'
echo ""

echo "IN CPP (first for double or int, second for string, third for char)"
echo '} else if (line.rfind("'$parname' =", 0) == 0 ){
'$parname'_ = std::stod(line.erase(0,line.find("=")+2));'
echo '} else if (line.rfind("'$parname' =", 0) == 0 ){
'$parname'_ = line.erase(0,line.find("=")+2);'
echo '} else if (line.rfind("'$parname' =", 0) == 0 ){
'$parname'_ = line.at(line.find("=")+2);'
echo 'stdoutfile << std::setw(8) << "'$parname':              "	<< '$parname'_ << std::endl;'
