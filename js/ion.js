function Ion (t, color, charge, terminus) {
	this.type = t;
	this.color = color;
	this.charge = charge;
	this.label = this.type;
	if(this.charge > 1)
		this.label += charge;
	this.label += "+";
	this.term = terminus;
}

// Source: http://en.wikipedia.org/wiki/Web_colors

// charge +1
Ion.A_1 = new Ion("a", "#008000", 1, "n"); // green
Ion.B_1 = new Ion("b", "#0000ff", 1, "n"); // blue
Ion.C_1 = new Ion("c", "#008B8B", 1, "n"); // dark cyan
Ion.X_1 = new Ion("x", "#4B0082", 1, "c"); // indigo
Ion.Y_1 = new Ion("y", "#ff0000", 1, "c"); // red
Ion.Z_1 = new Ion("z", "#FF8C00", 1, "c"); // dark orange

// charge +2
Ion.A_2 = new Ion("a", "#2E8B57", 2, "n"); // sea green
Ion.B_2 = new Ion("b", "#4169E1", 2, "n"); // royal blue
Ion.C_2 = new Ion("c", "#20B2AA", 2, "n"); // light sea green
Ion.X_2 = new Ion("x", "#800080", 2, "c"); // purple
Ion.Y_2 = new Ion("y", "#FA8072", 2, "c"); // salmon 
Ion.Z_2 = new Ion("z", "#FFA500", 2, "c"); // orange 

// charge +3
Ion.A_3 = new Ion("a", "#9ACD32", 3, "n"); // yellow green
Ion.B_3 = new Ion("b", "#00BFFF", 3, "n"); // deep sky blue
Ion.C_3 = new Ion("c", "#66CDAA", 3, "c"); // medium aquamarine
Ion.X_3 = new Ion("x", "#9932CC", 3, "c"); // dark orchid
Ion.Y_3 = new Ion("y", "#FFA07A", 3, "c"); // light salmon
Ion.Z_3 = new Ion("z", "#FFD700", 3, "n"); // gold

var _ions = [];
_ions["a"] = [];
_ions["a"][1] = Ion.A_1;
_ions["a"][2] = Ion.A_2;
_ions["a"][3] = Ion.A_3;
_ions["b"] = [];
_ions["b"][1] = Ion.B_1;
_ions["b"][2] = Ion.B_2;
_ions["b"][3] = Ion.B_3;
_ions["c"] = [];
_ions["c"][1] = Ion.C_1;
_ions["c"][2] = Ion.C_2;
_ions["c"][3] = Ion.C_3;
_ions["x"] = [];
_ions["x"][1] = Ion.X_1;
_ions["x"][2] = Ion.X_2;
_ions["x"][3] = Ion.X_3;
_ions["y"] = [];
_ions["y"][1] = Ion.Y_1;
_ions["y"][2] = Ion.Y_2;
_ions["y"][3] = Ion.Y_3;
_ions["z"] = [];
_ions["z"][1] = Ion.Z_1;
_ions["z"][2] = Ion.Z_2;
_ions["z"][3] = Ion.Z_3;

Ion.get = function _getIon(type, charge) {
	
	return _ions[type][charge];
}

Ion.getSeriesColor = function _getSeriesColor(ion) {
	
	return _ions[ion.type][ion.charge].color;
}


//-----------------------------------------------------------------------------
// Ion Series
//-----------------------------------------------------------------------------
var MASS_H = 1.00794;
var MASS_C = 12.011;
var MASS_N = 14.00674;
var MASS_O = 15.9994;
var MASS_PROTON = 1.007276;

Ion.MASS_PROTON = MASS_PROTON;
Ion.MASS_H = MASS_H;
Ion.MASS_C = MASS_C;
Ion.MASS_N = MASS_N;
Ion.MASS_O = MASS_O;

// massType can be "mono" or "avg"
Ion.getSeriesIon = function _getSeriesIon(ion, sequence, idxInSeq, massType) {
	if(ion.type == "a")	
		return new Ion_A (sequence, idxInSeq, ion.charge, massType);
	if(ion.type == "b")
		return new Ion_B (sequence, idxInSeq, ion.charge, massType);
	if(ion.type == "c")
		return new Ion_C (sequence, idxInSeq, ion.charge, massType);
	if(ion.type == "x")
		return new Ion_X (sequence, idxInSeq, ion.charge, massType);
	if(ion.type == "y")
		return new Ion_Y (sequence, idxInSeq, ion.charge, massType);
	if(ion.type == "z")
		return new Ion_Z (sequence, idxInSeq, ion.charge, massType);
}

function _makeIonLabel(type, index, charge) {
	var label = type+""+index;
	for(var i = 1; i <= charge; i+=1) 
		label += "+";
	return label;
}

function _getMz(neutralMass, charge) {
	return ( neutralMass + (charge * MASS_PROTON) ) / charge;
}

function _getWaterLossMz(sion) {
	var neutralMass = (sion.mz * sion.charge) - (sion.charge * MASS_PROTON);
	return _getMz((neutralMass - (MASS_H * 2 + MASS_O)), sion.charge);
}

function _getAmmoniaLossMz(sion) {
	var neutralMass = (sion.mz * sion.charge) - (sion.charge * MASS_PROTON);
	return _getMz((neutralMass - (MASS_H * 3 + MASS_N)), sion.charge);
}

Ion.getMz = _getMz;
Ion.getWaterLossMz = _getWaterLossMz;
Ion.getAmmoniaLossMz = _getAmmoniaLossMz;

function Ion_A (sequence, endIdxPlusOne, charge, massType) {
	// Neutral mass:  	 [N]+[M]-CHO  ; N = mass of neutral N terminal group
	var mass = 0;
	if(massType == "mono")
		mass = Peptide.getSeqMassMono(sequence, endIdxPlusOne, "n") - (MASS_C + MASS_O);
	else if(massType == "avg")
		mass = Peptide.getSeqMassAvg(sequence, endIdxPlusOne, "n") - (MASS_C + MASS_O);
	this.charge = charge;
	this.mz = _getMz(mass, charge);
	this.label = _makeIonLabel("a",endIdxPlusOne, charge);
	this.match = false;
	this.term = "n";
	return this;
}

function Ion_B (sequence, endIdxPlusOne, charge, massType) {
	// Neutral mass:    [N]+[M]-H  ; N = mass of neutral N terminal group
	var mass = 0;
	if(massType == "mono")
		mass = Peptide.getSeqMassMono(sequence, endIdxPlusOne, "n");
	else if(massType == "avg")
		mass = Peptide.getSeqMassAvg(sequence, endIdxPlusOne, "n");
	this.charge = charge;
	this.mz = _getMz(mass, charge);
	this.label = _makeIonLabel("b", endIdxPlusOne, charge);
	this.match = false;
	this.term = "n";
	return this;
}

function Ion_C (sequence, endIdxPlusOne, charge, massType) {
	// Neutral mass:    [N]+[M]+NH2  ; N = mass of neutral N terminal group
	var mass = 0;
	if(massType == "mono")
		mass = Peptide.getSeqMassMono(sequence, endIdxPlusOne, "n") + MASS_H + (MASS_N + 2*MASS_H);
	else if(massType == "avg")
		mass = Peptide.getSeqMassAvg(sequence, endIdxPlusOne, "n") + MASS_H + (MASS_N + 2*MASS_H);
	this.charge = charge;
	this.mz = _getMz(mass, charge);
	this.label = _makeIonLabel("c", endIdxPlusOne, charge);
	this.match = false;
	this.term = "n";
	return this;
}

function Ion_X (sequence, startIdx, charge, massType) {
	// Neutral mass = [C]+[M]+CO-H ; C = mass of neutral C-terminal group (OH)
	var mass = 0;
	if(massType == "mono")
		mass = Peptide.getSeqMassMono(sequence, startIdx, "c") + 2*MASS_O + MASS_C;
	else if(massType == "avg")
		mass = Peptide.getSeqMassAvg(sequence, startIdx, "c") + 2*MASS_O + MASS_C;
	this.charge = charge;
	this.mz = _getMz(mass, charge);
	this.label = _makeIonLabel("x", sequence.length - startIdx, charge);
	this.match = false;
	this.term = "c";
	return this;
}

function Ion_Y (sequence, startIdx, charge, massType) {
	// Neutral mass = [C]+[M]+H ; C = mass of neutral C-terminal group (OH)
	var mass = 0;
	if(massType == "mono")
		mass = Peptide.getSeqMassMono(sequence, startIdx, "c") + 2*MASS_H + MASS_O;
	else if(massType == "avg")
		mass = Peptide.getSeqMassAvg(sequence, startIdx, "c") + 2*MASS_H + MASS_O;
	this.charge = charge;
	this.mz = _getMz(mass, charge);
	this.label = _makeIonLabel("y", sequence.length - startIdx, charge);
	this.match = false;
	this.term = "c";
	return this;
}

function Ion_Z (sequence, startIdx, charge, massType) {
	// Neutral mass = [C]+[M]-NH2 ; C = mass of neutral C-terminal group (OH)
	// We're really printing Z-dot ions so we add an H to make it OH+[M]-NH2 +H = [M]+O-N
	var mass = 0;
	if(massType == "mono")
		mass = Peptide.getSeqMassMono(sequence, startIdx, "c") + MASS_O - MASS_N;
	else if(massType == "avg")
		mass = Peptide.getSeqMassAvg(sequence, startIdx, "c") + MASS_O - MASS_N;
	this.charge = charge;
	this.mz = _getMz(mass, charge);
	this.label = _makeIonLabel("z", sequence.length - startIdx, charge);
	this.match = false;
	this.term = "c";
	return this;
}