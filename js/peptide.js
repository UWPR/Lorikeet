// -----------------------------------------------------------------------------
// Peptide sequence and modifications
// -----------------------------------------------------------------------------
Peptide.sequence;
Peptide.staticMods;
Peptide.varMods;
Peptide.ntermMod;
Peptide.ctermMod;

function Peptide(seq, staticModifications, varModifications, ntermModification, ctermModification) {
	
	Peptide.sequence = seq;
	Peptide.ntermMod = ntermModification;
	Peptide.ctermMod = ctermModification;
	Peptide.staticMods = [];
	if(staticModifications) {
		for(var i = 0; i < staticModifications.length; i += 1) {
			var mod = staticModifications[i];
			Peptide.staticMods[mod.aa.code] = mod;
		}
	}
	
	Peptide.varMods = [];
	if(varModifications) {
		for(var i = 0; i < varModifications.length; i += 1) {
			var mod = varModifications[i];
			Peptide.varMods[mod.position] = mod;
		}
	}
}


// index: index in the seq.
// If this is a N-term sequence we will sum up the mass of the amino acids in the sequence up-to index (exclusive).
// If this is a C-term sequence we will sum up the mass of the amino acids in the sequence starting from index (inclusive)
// modification masses are added
Peptide.getSeqMassMono = function _seqMassMono(seq, index, term) {
	
	var mass = 0;
	var aa_obj = new AminoAcid();
	if(seq) {
		if(term == "n") {
			for( var i = 0; i < index; i += 1) {
				var aa = aa_obj.get(seq[i]);
				mass += aa.mono;
			}
		}
		if (term == "c") {
			for( var i = index; i < seq.length; i += 1) {
				var aa = aa_obj.get(seq[i]);
				mass += aa.mono;
			}
		}
	}
	
	mass = _addTerminalModMass(mass, term);
	mass = _addResidueModMasses(mass, seq, index, term);
	return mass;
}

// index: index in the seq.
// If this is a N-term sequence we will sum up the mass of the amino acids in the sequence up-to index (exclusive).
// If this is a C-term sequence we will sum up the mass of the amino acids in the sequence starting from index (inclusive)
// modification masses are added
Peptide.getSeqMassAvg = function _seqMassAvg(seq, index, term) {
	
	var mass = 0;
	var aa_obj = new AminoAcid();
	if(seq) {
		if(term == "n") {
			for( var i = 0; i < index; i += 1) {
				var aa = aa_obj.get(seq[i]);
				mass += aa.avg;
			}
		}
		if (term == "c") {
			for( var i = index; i < seq.length; i += 1) {
				var aa = aa_obj.get(seq[i]);
				mass += aa.avg;
			}
		}
	}
	
	mass = _addTerminalModMass(mass, term);
	mass = _addResidueModMasses(mass, seq, index, term);
	return mass;
}

// Returns the monoisotopic neutral mass of the peptide; modifications added. N-term H and C-term OH are added
Peptide.getNeutralMassMono = function _massNeutralMono() {
	
	var mass = 0;
	var aa_obj = new AminoAcid();
	if(Peptide.sequence) {
		for(var i = 0; i < Peptide.sequence.length; i++) {
			var aa = aa_obj.get(Peptide.sequence[i]);
			mass += aa.mono;
		}
	}
	
	mass = _addTerminalModMass(mass, "n");
	mass = _addTerminalModMass(mass, "c");
	mass = _addResidueModMasses(mass, Peptide.sequence, Peptide.sequence.length, "n");
	// add N-terminal H
	mass = mass + Ion.MASS_H_1;
	// add C-terminal OH
	mass = mass + Ion.MASS_O_16 + Ion.MASS_H_1;
	
	return mass;
}

//Returns the avg neutral mass of the peptide; modifications added. N-term H and C-term OH are added
Peptide.getNeutralMassAvg = function _massNeutralAvg() {
	
	var mass = 0;
	var aa_obj = new AminoAcid();
	if(Peptide.sequence) {
		for(var i = 0; i < Peptide.sequence.length; i++) {
			var aa = aa_obj.get(Peptide.sequence[i]);
			mass += aa.avg;
		}
	}
	
	mass = _addTerminalModMass(mass, "n");
	mass = _addTerminalModMass(mass, "c");
	mass = _addResidueModMasses(mass, Peptide.sequence, Peptide.sequence.length, "n");
	// add N-terminal H
	mass = mass + Ion.MASS_H;
	// add C-terminal OH
	mass = mass + Ion.MASS_O + Ion.MASS_H;
	
	return mass;
}


function _addResidueModMasses(seqMass, seq, index, term) {
	
	var mass = seqMass;
	
	// add any static modifications
	if(term == "n") {
		for(var i = 0; i < index; i += 1) {
			var mod = Peptide.staticMods[seq[i]];
			if(mod) {
				mass += mod.modMass;
			}
		}
	}
	if(term == "c") {
		for(var i = index; i < seq.length; i += 1) {
			var mod = Peptide.staticMods[seq[i]];
			if(mod) {
				mass += mod.modMass;
			}
		}
	}
	
	// add any variable modifications
	if(term == "n") {
		for(var i = 0; i < index; i += 1) {
			var mod = Peptide.varMods[i+1]; // varMods index in the sequence is 1-based
			if(mod) {
				mass += mod.modMass;
			}
		}
	}
	if(term == "c") {
		for(var i = index; i < seq.length; i += 1) {
			var mod = Peptide.varMods[i+1]; // varMods index in the sequence is 1-based
			if(mod) {
				mass += mod.modMass;
			}
		}
	}
	return mass;
}

function _addTerminalModMass(seqMass, term) {
	
	var mass = seqMass;
	// add any terminal modifications
	if(term == "n" && Peptide.ntermMod)
		mass += Peptide.ntermMod;
	if(term == "c" && Peptide.ctermMod)
		mass += Peptide.ctermMod;

	return mass;
}


//-----------------------------------------------------------------------------
// Modification
//-----------------------------------------------------------------------------
function Modification(aminoAcid, mass) {
	this.aa = aminoAcid;
	this.modMass = mass;
}

function VariableModification(pos, mass, aminoAcid) {
	this.position = parseInt(pos);
	this.aa = aminoAcid;
	this.modMass = mass;
}



