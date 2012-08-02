echo "
/* File: smac.i */
%module smac
%{
	#include \"../../src/smac.h\"
%}

%include <std_vector.i>
%include \"std_complex.i\"

namespace std
{
	%template(int_array_t) vector<int>;
	%template(bool_array_t) vector<bool>;
	%template(coord_t) vector<double>;
	%template(coordlist_t) vector< vector<double> >;
	%template(component_t) complex<double>;	
	%template(shpdesc_t) vector <complex<double> >;
}
"

for i in $(ls ../../src/*.h); do
	echo "%include \"$i\""
done
