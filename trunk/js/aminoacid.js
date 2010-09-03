
function AminoAcid(aaCode, aaShortName, aaName, monoMass, avgMass) {
	this.code = aaCode;
	this.shortName = aaShortName;
	this.name = aaName;
	this.mono = monoMass;
	this.avg = avgMass;
	
	this.get = _getAA;
}

AminoAcid.A = new AminoAcid ("A", "Ala", "Alanine", 71.03711, 71.0788);
AminoAcid.R = new AminoAcid ("R", "Arg", "Arginine", 156.10111, 156.1875);
AminoAcid.N = new AminoAcid ("N", "Asn", "Asparagine", 114.04293, 114.1038);
AminoAcid.D = new AminoAcid ("D", "Asp", "Aspartic Acid", 115.02694, 115.0886);
AminoAcid.C = new AminoAcid ("C", "Cys", "Cysteine", 103.00919, 103.1388);
AminoAcid.E = new AminoAcid ("E", "Glu", "Glutamine", 129.04259, 129.1155);
AminoAcid.Q = new AminoAcid ("Q", "Gln", "Glutamic Acid", 128.05858, 128.1307);
AminoAcid.G = new AminoAcid ("G", "Gly", "Glycine", 57.02146, 57.0519);
AminoAcid.H = new AminoAcid ("H", "His", "Histidine", 137.05891, 137.1411);
AminoAcid.I = new AminoAcid ("I", "Ile", "Isoleucine", 113.08406, 113.1594);
AminoAcid.L = new AminoAcid ("L", "Leu", "Leucine", 113.08406, 113.1594);
AminoAcid.K = new AminoAcid ("K", "Lys", "Lysine", 128.09496, 128.1741);
AminoAcid.M = new AminoAcid ("M", "Met", "Methionine", 131.04049, 131.1926);
AminoAcid.F = new AminoAcid ("F", "Phe", "Phenylalanine", 147.06841, 147.1766);
AminoAcid.P = new AminoAcid ("P", "Pro", "Proline", 97.05276, 97.1167);
AminoAcid.S = new AminoAcid ("S", "Ser", "Serine", 87.03203, 87.0782);
AminoAcid.T = new AminoAcid ("T", "Thr", "Threonine", 101.04768, 101.1051);
AminoAcid.W = new AminoAcid ("W", "Trp", "Tryptophan", 186.07931, 186.2132);
AminoAcid.Y = new AminoAcid ("Y", "Tyr", "Tyrosine", 163.06333, 163.1760);
AminoAcid.V = new AminoAcid ("V", "Val", "Valine", 99.06841, 99.1326);


AminoAcid.aa = [];
AminoAcid.aa["A"] = AminoAcid.A;
AminoAcid.aa["R"] = AminoAcid.R;
AminoAcid.aa["N"] = AminoAcid.N;
AminoAcid.aa["D"] = AminoAcid.D;
AminoAcid.aa["C"] = AminoAcid.C;
AminoAcid.aa["E"] = AminoAcid.E;
AminoAcid.aa["Q"] = AminoAcid.Q;
AminoAcid.aa["G"] = AminoAcid.G;
AminoAcid.aa["H"] = AminoAcid.H;
AminoAcid.aa["I"] = AminoAcid.I;
AminoAcid.aa["L"] = AminoAcid.L;
AminoAcid.aa["K"] = AminoAcid.K;
AminoAcid.aa["M"] = AminoAcid.M;
AminoAcid.aa["F"] = AminoAcid.F;
AminoAcid.aa["P"] = AminoAcid.P;
AminoAcid.aa["S"] = AminoAcid.S;
AminoAcid.aa["T"] = AminoAcid.T;
AminoAcid.aa["W"] = AminoAcid.W;
AminoAcid.aa["Y"] = AminoAcid.Y;
AminoAcid.aa["V"] = AminoAcid.V;

AminoAcid.get = _getAA;

function _getAA(aaCode) {
	if(AminoAcid.aa[aaCode])
		return AminoAcid.aa[aaCode];
	else
		return new AminoAcid(aaCode, aaCode, 0.0, 0.0);
}
