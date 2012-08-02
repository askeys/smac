
/* File: smac.i */
%module smac
%{
	#include "../../src/smac.h"
	using namespace std;
	using namespace smac;
%}

%include <std_vector.i>
%include <std_complex.i>
%include <std_string.i>
%include <carrays.i>
%include <cdata.i>

namespace std
{
	%template(int_array) vector<int>;
	%template(bool_array) vector<bool>;
	%template(coord) vector<double>;
	%template(coordlist) vector< vector<double> >;
//	%template(component) std::complex<double>;
	%template(string_array) vector<string>;
	%template(shpdesc) vector <complex<double> >;
}

%include "../../src/Cell.h"
%include "../../src/CellList.h"
%include "../../src/cluster.h"
%include "../../src/config.h"
%include "../../src/diagnostics.h"
%include "../../src/io.h"
%include "../../src/matchfun.h"
%include "../../src/neighbormap.h"
%include "../../src/operator.h"
%include "../../src/preprocess.h"
%include "../../src/rand.h"
%include "../../src/register.h"
%include "../../src/shpdesccontext.h"
%include "../../src/shpdescdist.h"
%include "../../src/shpdescfourier.h"
%include "../../src/shpdeschist.h"
%include "../../src/shpdescrms.h"
%include "../../src/shpdescwigner3j.h"
%include "../../src/shpdesczernike.h"
%include "../../src/smac.h"
%include "../../src/space.h"
%include "../../src/stdanalysis.h"
%include "../../src/typedefs.h"
%include "../../src/voxel.h"
