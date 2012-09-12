(function($) {

    // plugin name - specview
	$.fn.specview = function (opts) {

        var defaults = {
                sequence: null,
                scanNum: null,
                fileName: null,
                charge: null,
                precursorMz: null,
                staticMods: [],
                variableMods: [],
                ntermMod: 0, // additional mass to be added to the n-term
                ctermMod: 0, // additional mass to be added to the c-term
                peaks: [],
                sparsePeaks: null,
                ms1peaks: null,
                ms1scanLabel: null,
                precursorPeaks: null,
                precursorPeakClickFn: null,
                zoomMs1: false,
                width: 750, 	// width of the ms/ms plot
                height: 450, 	// height of the ms/ms plot
                massError: 0.5, // mass tolerance for labeling peaks
                extraPeakSeries:[],
                showIonTable: true,
                showViewingOptions: true
        };
			
	    var options = $.extend(true, {}, defaults, opts); // this is a deep copy

        return this.each(function() {

            index = index + 1;
            init($(this), options);

        });
	};

    var index = 0;

    var elementIds = {
            massError: "massError",
            msPlot: "msPlot",
            msmsplot: "msmsplot",
            ms1plot_zoom_out: "ms1plot_zoom_out",
            ms1plot_zoom_in: "ms1plot_zoom_in",
            ms2plot_zoom_out: "ms2plot_zoom_out",
            zoom_x: "zoom_x",
            zoom_y: "zoom_y",
            resetZoom: "resetZoom",
            update: "update",
            enableTooltip: "enableTooltip",
            msmstooltip: "lorikeet_msmstooltip",
            ion_choice: "ion_choice",
            nl_choice: "nl_choice",
            deselectIonsLink: "deselectIonsLink",
            slider_width: "slider_width",
            slider_width_val: "slider_width_val",
            slider_height: "slider_height",
            slider_height_val: "slider_height_val",
            printLink: "printLink",
            lorikeet_content: "lorikeet_content",
            optionsTable: "optionsTable",
            ionTableLoc1: "ionTableLoc1",
            ionTableLoc2: "ionTableLoc2",
            viewOptionsDiv: "viewOptionsDiv",
            moveIonTable: "moveIonTable",
            modInfo: "modInfo",
            ionTableDiv: "ionTableDiv",
            ionTable: "ionTable",
            fileinfo: "fileinfo",
            seqinfo: "seqinfo",
            removeNoise: "removeNoise"
	};

    function getElementId(container, elementId){
        return elementId+"_"+container.data("index");
    }

    function getElementSelector(container, elementId) {
        return "#"+getElementId(container, elementId);
    }

    function getRadioName(container, name) {
        return name+"_"+container.data("index");
    }

    function init(parent_container, options) {

        // read the static modifications
        var parsedStaticMods = [];
        for(var i = 0; i < options.staticMods.length; i += 1) {
            var mod = options.staticMods[i];
            parsedStaticMods[i] = new Modification(AminoAcid.get(mod.aminoAcid), mod.modMass);
        }
        options.staticMods = parsedStaticMods;

        // read the variable modifications
        var parsedVarMods = [];
        for(var i = 0; i < options.variableMods.length; i += 1) {
            // position: 14, modMass: 16.0, aminoAcid: 'M'
            var mod = options.variableMods[i];
            parsedVarMods[i] = new VariableModification(
                                    mod.index,
                                    mod.modMass,
                                    AminoAcid.get(mod.aminoAcid)
                                );
        }
        options.variableMods = parsedVarMods;

        var peptide = new Peptide(options.sequence, options.staticMods, options.variableMods,
                                options.ntermMod, options.ctermMod);
        options.peptide = peptide;

        var container = createContainer(parent_container);
        // alert(container.attr('id')+" parent "+container.parent().attr('id'));
        storeContainerData(container, options);
        initContainer(container);

        if(options.charge) {
            makeOptionsTable(container, Math.max(1,options.charge-1));
        }
        else
            makeOptionsTable(container, 1);

        makeViewingOptions(container, options);

        showSequenceInfo(container, options);
        showFileInfo(container, options);
        showModInfo(container, options);


        createPlot(container, getDatasets(container)); // Initial MS/MS Plot
        if(options.ms1peaks && options.ms1peaks.length > 0) {
            if(options.zoomMs1 && options.precursorMz) {

                var ms1zoomRange = container.data("ms1zoomRange");

                ms1zoomRange = {xaxis: {}, yaxis: {}};
                ms1zoomRange.xaxis.from = options.precursorMz - 5.0;
                ms1zoomRange.xaxis.to = options.precursorMz + 5.0;
                var max_intensity = 0;
                for(var j = 0; j < options.ms1peaks.length; j += 1) {
                    var pk = options.ms1peaks[j];
                    if(pk[0] < options.precursorMz - 5.0)
                        continue;
                    if(pk[0] > options.precursorMz + 5.0)
                        break;
                    if(pk[1] > max_intensity)
                        max_intensity = pk[1];
                }
                ms1zoomRange.yaxis.from = 0.0;
                ms1zoomRange.yaxis.to = max_intensity;
            }
            createMs1Plot(container);
            setupMs1PlotInteractions(container);
        }

        setupInteractions(container);

        makeIonTable(container);
    }

    function storeContainerData(container, options) {

        container.data("index", index);
        container.data("options", options);
        container.data("massErrorChanged", false);
        container.data("massTypeChanged", false);
        container.data("peakAssignmentTypeChanged", false);
        container.data("peakLabelTypeChanged", false);
        container.data("selectedNeutralLossChanged", false);
        container.data("plot", null);           // MS/MS plot
        container.data("ms1plot", null);        // MS1 plot (only created when data is available)
        container.data("zoomRange", null);      // for zooming MS/MS plot
        container.data("ms1zoomRange", null);
        container.data("previousPoint", null);  // for tooltips
        container.data("ionSeries", {a: [], b: [], c: [], x: [], y: [], z: []});
        container.data("ionSeriesLabels", {a: [], b: [], c: [], x: [], y: [], z: []});
        container.data("ionSeriesMatch", {a: [], b: [], c: [], x: [], y: [], z: []});
        container.data("massError", options.massError);

        var maxInt = getMaxInt(options);
        var plotOptions =  {
                series: {
                    peaks: { show: true, lineWidth: 1, shadowSize: 0},
                    shadowSize: 0
                },
                selection: { mode: "x", color: "#F0E68C" },
                grid: { show: true,
                        hoverable: true,
                        clickable: false,
                        autoHighlight: false,
                        borderWidth: 1,
                        labelMargin: 1},
                xaxis: { tickLength: 3, tickColor: "#000" },
                yaxis: { tickLength: 0, tickColor: "#000",
                        tickFormatter: function(val, axis) {return Math.round((val * 100)/maxInt)+"%";}}
	        }
        container.data("plotOptions", plotOptions);
        container.data("maxInt", maxInt);

    }

	function getMaxInt(options) {
		var maxInt = 0;
		for(var j = 0; j < options.peaks.length; j += 1) {
			var peak = options.peaks[j];
			if(peak[1] > maxInt) {
				maxInt = peak[1];
			}
		}
		//alert(maxInt);
		return maxInt;
	}
	
	function round(number) {
		return number.toFixed(4);
	}
	
	function setMassError(container) {
   		var me = parseFloat($(getElementSelector(container, elementIds.massError)).val());
		if(me != container.data("massError")) {
			container.data("massError", me);
			container.data("massErrorChanged", true);
		}
		else {
			container.data("massErrorChanged", false);
		}
   	}

    // -----------------------------------------------
	// CREATE MS1 PLOT
	// -----------------------------------------------
	function createMs1Plot(container) {

        var ms1zoomRange = container.data("ms1zoomRange");
        var options = container.data("options");

		var data = [{data: options.ms1peaks, color: "#bbbbbb", labelType: 'none', hoverable: false, clickable: false}];
		if(options.precursorPeaks) {
			if(options.precursorPeakClickFn)
				data.push({data: options.precursorPeaks, color: "#ff0000", hoverable: true, clickable: true});
			else
				data.push({data: options.precursorPeaks, color: "#ff0000", hoverable: false, clickable: false});
		}
		
		// the MS/MS plot should have been created by now.  This is a hack to get the plots aligned.
		// We will set the y-axis labelWidth to this value.
		var labelWidth = container.data("plot").getAxes().yaxis.labelWidth;
		
		var ms1plotOptions = {
				series: { peaks: {show: true, shadowSize: 0}, shadowSize: 0},
				grid: { show: true, 
						hoverable: true, 
						autoHighlight: true,
						clickable: true,
						borderWidth: 1,
						labelMargin: 1},
				selection: { mode: "xy", color: "#F0E68C" },
		        xaxis: { tickLength: 2, tickColor: "#000" },
		    	yaxis: { tickLength: 0, tickColor: "#fff", ticks: [], labelWidth: labelWidth }
		};
		
		if(ms1zoomRange) {
			ms1plotOptions.xaxis.min = ms1zoomRange.xaxis.from;
			ms1plotOptions.xaxis.max = ms1zoomRange.xaxis.to;
			ms1plotOptions.yaxis.min = ms1zoomRange.yaxis.from;
			ms1plotOptions.yaxis.max = ms1zoomRange.yaxis.to;
		}

		var placeholder = $(getElementSelector(container, elementIds.msPlot));
		var ms1plot = $.plot(placeholder, data, ms1plotOptions);
        container.data("ms1plot", ms1plot);


		// mark the current precursor peak
		if(options.precursorPeaks) {
			var x,y, diff, precursorMz;
			
			// If we are given a precursor m/z use it
			if(options.precursorMz) {
				precursorMz = options.precursorMz;
			}
			// Otherwise calculate a theoretical m/z from the given sequence and charge
			else if(options.sequence && options.charge) {
				var mass = options.peptide.getNeutralMassMono();
				precursorMz = Ion.getMz(mass, options.charge);
			}
			
			if(precursorMz) {
				// find the closest actual peak
				for(var i = 0; i < options.precursorPeaks.length; i += 1) {
					var pk = options.precursorPeaks[i];
					var d = Math.abs(pk[0] - precursorMz);
					if(!diff || d < diff) {
						x = pk[0];
						y = pk[1];
						diff = d;
					}
				}
				if(diff <= 0.5) {
					var o = ms1plot.pointOffset({ x: x, y: y});
				    var ctx = ms1plot.getCanvas().getContext("2d");
				    ctx.beginPath();
				    ctx.moveTo(o.left-10, o.top-5);
				    ctx.lineTo(o.left-10, o.top + 5);
				    ctx.lineTo(o.left-10 + 10, o.top);
				    ctx.lineTo(o.left-10, o.top-5);
				    ctx.fillStyle = "#008800";
				    ctx.fill();
				    placeholder.append('<div style="position:absolute;left:' + (o.left + 4) + 'px;top:' + (o.top-4) + 'px;color:#000;font-size:smaller">'+x.toFixed(2)+'</div>');
				}
			}
		}
		
		// mark the scan number if we have it
		o = ms1plot.getPlotOffset();
		if(options.ms1scanLabel) {
			placeholder.append('<div style="position:absolute;left:' + (o.left + 4) + 'px;top:' + (o.top+4) + 'px;color:#666;font-size:smaller">MS1 scan: '+options.ms1scanLabel+'</div>');
		}
		
		// zoom out icon on plot right hand corner
		if(container.data("zoomrange")) {
			placeholder.append('<div id="'+getElementId(container, elementIds.ms1plot_zoom_out)+'" class="zoom_out_link"  style="position:absolute; left:'
					+ (o.left + ms1plot.width() - 40) + 'px;top:' + (o.top+4) + 'px;"></div>');

			$(getElementSelector(container, elementIds.ms1plot_zoom_out)).click( function() {
				ms1zoomRange = null;
				createMs1Plot(container);
			});
		}
		
		if(options.precursorPeaks) {
			placeholder.append('<div id="'+getElementId(container, elementIds.ms1plot_zoom_in)+'" class="zoom_in_link"  style="position:absolute; left:'
					+ (o.left + ms1plot.width() - 20) + 'px;top:' + (o.top+4) + 'px;"></div>');
			$(getElementSelector(container, elementIds.ms1plot_zoom_in)).click( function() {
				var ranges = {};
				ranges.yaxis = {};
				ranges.xaxis = {};
				ranges.yaxis.from = null;
				ranges.yaxis.to = null;
				ranges.xaxis.from = null;
				ranges.xaxis.to = null;
				var maxInt = 0;
				for(var p = 0; p < options.precursorPeaks.length; p += 1) {
					if(options.precursorPeaks[p][1] > maxInt)
						maxInt = options.precursorPeaks[p][1];
				}
				ranges.yaxis.to = maxInt;
                container.data("ms1zoomRange", ranges);
				createMs1Plot(container);
			});
		}
		

	}

    // -----------------------------------------------
	// SET UP INTERACTIVE ACTIONS FOR MS1 PLOT
	// -----------------------------------------------
	function setupMs1PlotInteractions(container) {
		
		var placeholder = $(getElementSelector(container, elementIds.msPlot));
        var options = container.data("options");

		// allow clicking on plot if we have a function to handle the click
		if(options.precursorPeakClickFn != null) {
			placeholder.bind("plotclick", function (event, pos, item) {
				
		        if (item) {
		          //highlight(item.series, item.datapoint);
		        	options.precursorPeakClickFn(item.datapoint[0]);
		        }
		    });
		}
		
		// allow zooming the plot
		placeholder.bind("plotselected", function (event, ranges) {
            container.data("ms1zoomRange", ranges);
			createMs1Plot(container);
	    });
		
	}

    // -----------------------------------------------
	// CREATE MS/MS PLOT
	// -----------------------------------------------
	function createPlot(container, datasets) {

        var plot;
    	if(!container.data("zoomRange")) {
    		plot = $.plot($(getElementSelector(container, elementIds.msmsplot)), datasets,  container.data("plotOptions"));
        }
    	else {
            var zoomRange = container.data("zoomRange");
            var selectOpts = {};
    		if($(getElementSelector(container, elementIds.zoom_x)).is(":checked"))
    			selectOpts.xaxis = { min: zoomRange.xaxis.from, max: zoomRange.xaxis.to };
    		if($(getElementSelector(container, elementIds.zoom_y)).is(":checked"))
    			selectOpts.yaxis = { min: zoomRange.yaxis.from, max: zoomRange.yaxis.to };
    		
    		plot = $.plot(getElementSelector(container, elementIds.msmsplot), datasets,
                      $.extend(true, {}, container.data("plotOptions"), selectOpts));

    		// zoom out icon on plot right hand corner
    		var o = plot.getPlotOffset();
    		$(getElementSelector(container, elementIds.msmsplot)).append('<div id="'+getElementId(container, elementIds.ms2plot_zoom_out)+'" class="zoom_out_link" style="position:absolute; left:'
					+ (o.left + plot.width() - 20) + 'px;top:' + (o.top+4) + 'px"></div>');

			$(getElementSelector(container, elementIds.ms2plot_zoom_out)).click( function() {
                resetZoom(container);
			});
    	}
    	// we have re-calculated and re-drawn everything..
    	container.data("massTypeChanged", false);
    	container.data("massErrorChanged",false);
    	container.data("peakAssignmentTypeChanged", false);
    	container.data("peakLabelTypeChanged", false);
    	container.data("selectedNeutralLossChanged", false);
        container.data("plot", plot);
    }
	
	// -----------------------------------------------
	// SET UP INTERACTIVE ACTIONS FOR MS/MS PLOT
	// -----------------------------------------------
	function setupInteractions (container) {

		// ZOOMING
	    $(getElementSelector(container, elementIds.msmsplot)).bind("plotselected", function (event, ranges) {
	    	container.data("zoomRange", ranges);
	    	createPlot(container, getDatasets(container));
	    });
	    
	    // ZOOM AXES
	    $(getElementSelector(container, elementIds.zoom_x)).click(function() {
	    	resetAxisZoom(container);
	    });
	    $(getElementSelector(container, elementIds.zoom_y)).click(function() {
	    	resetAxisZoom(container);
	    });
	    
		// RESET ZOOM
		$(getElementSelector(container, elementIds.resetZoom)).click(function() {
			resetZoom(container);
	   	});
		
		// UPDATE
		$(getElementSelector(container, elementIds.update)).click(function() {
			container.data("zoomRange", null); // zoom out fully
			setMassError(container);
			createPlot(container, getDatasets(container));
			makeIonTable(container);
	   	});
		
		// TOOLTIPS
		$(getElementSelector(container, elementIds.msmsplot)).bind("plothover", function (event, pos, item) {

	        if($(getElementSelector(container, elementIds.enableTooltip)+":checked").length > 0) {
	            if (item) {
	                if (container.data("previousPoint") != item.datapoint) {
	                    container.data("previousPoint", item.datapoint);
	                    
	                    $(getElementSelector(container, elementIds.msmstooltip)).remove();
	                    var x = item.datapoint[0].toFixed(2),
	                        y = item.datapoint[1].toFixed(2);
	                    
	                    showTooltip(container, item.pageX, item.pageY,
	                                "m/z: " + x + "<br>intensity: " + y);
	                }
	            }
	            else {
	                $(getElementSelector(container, elementIds.msmstooltip)).remove();
	                container.data("previousPoint", null);
	            }
	        }
	    });
		$(getElementSelector(container, elementIds.enableTooltip)).click(function() {
			$(getElementSelector(container, elementIds.msmstooltip)).remove();
		});
		
		// SHOW / HIDE ION SERIES; UPDATE ON MASS TYPE CHANGE; 
		// PEAK ASSIGNMENT TYPE CHANGED; PEAK LABEL TYPE CHANGED
		var ionChoiceContainer = $(getElementSelector(container, elementIds.ion_choice));
		ionChoiceContainer.find("input").click(function() {
            plotAccordingToChoices(container)
        });
		
		var neutralLossContainer = $(getElementSelector(container, elementIds.nl_choice));
		neutralLossContainer.find("input").click(function() {
			container.data("selectedNeutralLossChanged", true);
			plotAccordingToChoices(container);
		});
		
	    container.find("input[name='"+getRadioName(container, "massTypeOpt")+"']").click(function() {
	    	container.data("massTypeChanged", true);
	    	plotAccordingToChoices(container);
	    });
        $(getElementSelector(container, elementIds.removeNoise)).click(function() {
            container.data("peakAssignmentTypeChanged", true);
            plotAccordingToChoices(container);
        });

	    container.find("input[name='"+getRadioName(container, "peakAssignOpt")+"']").click(function() {
	    	container.data("peakAssignmentTypeChanged", true);
	    	plotAccordingToChoices(container);
	    });

        $(getElementSelector(container, elementIds.deselectIonsLink)).click(function() {
			ionChoiceContainer.find("input:checkbox:checked").each(function() {
				$(this).attr('checked', "");
			});
			
			plotAccordingToChoices(container);
		});

	    container.find("input[name='"+getRadioName(container, "peakLabelOpt")+"']").click(function() {
	    	container.data("peakLabelTypeChanged", true);
	    	plotAccordingToChoices(container);
	    });

	    
	    
	    // MOVING THE ION TABLE
	    makeIonTableMovable(container);
	    
	    // CHANGING THE PLOT SIZE
	    makePlotResizable(container);
	    
	    // PRINT SPECTRUM
	    printPlot(container);
		
	}
	
	function resetZoom(container) {
		container.data("zoomRange", null);
		setMassError(container);
		createPlot(container, getDatasets(container));
	}
	
	function plotAccordingToChoices(container) {
        var data = getDatasets(container);

		if (data.length > 0) {
            createPlot(container, data);
            makeIonTable(container);
            showSequenceInfo(container); // update the MH+ and m/z values
        }
    }
	
	function resetAxisZoom(container) {

        var plot = container.data("plot");
        var plotOptions = container.data("plotOptions");

    	var zoom_x = false;
		var zoom_y = false;
		if($(getElementSelector(container, elementIds.zoom_x)).is(":checked"))
			zoom_x = true;
		if($(getElementSelector(container, elementIds.zoom_y)).is(":checked"))
			zoom_y = true;
    	if(zoom_x && zoom_y) {
    		plotOptions.selection.mode = "xy";
			if(plot) plot.getOptions().selection.mode = "xy";
    	}
		else if(zoom_x) {
			plotOptions.selection.mode = "x";
			if(plot) plot.getOptions().selection.mode = "x";
		}
		else if(zoom_y) {
			plotOptions.selection.mode = "y";
			if(plot) plot.getOptions().selection.mode = "y";
		}
	}
	
	function showTooltip(container, x, y, contents) {
        $('<div id="'+getElementId(container, elementIds.msmstooltip)+'">' + contents + '</div>').css( {
            position: 'absolute',
            display: 'none',
            top: y + 5,
            left: x + 5,
            border: '1px solid #fdd',
            padding: '2px',
            'background-color': '#F0E68C',
            opacity: 0.80
        }).appendTo("body").fadeIn(200);
    }
	
	function makePlotResizable(container) {

        var options = container.data("options");

		$(getElementSelector(container, elementIds.slider_width)).slider({
			value:options.width,
			min: 100,
			max: 1500,
			step: 50,
			slide: function(event, ui) {
				var width = ui.value;
				//console.log(ui.value);
				options.width = width;
				$(getElementSelector(container, elementIds.msmsplot)).css({width: width});
				plotAccordingToChoices(container);
				if(options.ms1peaks && options.ms1peaks.length > 0) {
					$(getElementSelector(container, elementIds.msPlot)).css({width: width});
					createMs1Plot(ms1zoomRange);
				}
				$(getElementSelector(container, elementIds.slider_width_val)).text(width);
			}
		});
		
		$(getElementSelector(container, elementIds.slider_height)).slider({
			value:options.height,
			min: 100,
			max: 1000,
			step: 50,
			slide: function(event, ui) {
				var height = ui.value;
				//console.log(ui.value);
				options.height = height
				$(getElementSelector(container, elementIds.msmsplot)).css({height: height});
				plotAccordingToChoices(container);
				$(getElementSelector(container, elementIds.slider_height_val)).text(height);
			}
		});
	}
	
	function printPlot(container) {

		$(getElementSelector(container, elementIds.printLink)).click(function() {

            var parent = container.parent();

			// create another div and move the plots into that div
			$(document.body).append('<div id="tempPrintDiv"></div>');
			$("#tempPrintDiv").append(container.detach());
			$("#tempPrintDiv").siblings().addClass("noprint");
			
			var plotOptions = container.data("plotOptions");

			container.find(".bar").addClass('noprint');
			$(getElementSelector(container, elementIds.optionsTable)).addClass('noprint');
			$(getElementSelector(container, elementIds.ionTableLoc1)).addClass('noprint');
			$(getElementSelector(container, elementIds.ionTableLoc2)).addClass('noprint');
			$(getElementSelector(container, elementIds.viewOptionsDiv)).addClass('noprint');
			
			plotOptions.series.peaks.print = true; // draw the labels in the DOM for sharper print output
			plotAccordingToChoices(container);
			window.print();
			
			
			// remove the class after printing so that if the user prints 
			// via the browser's print menu the whole page is printed
			container.find(".bar").removeClass('noprint');
			$(getElementSelector(container, elementIds.optionsTable)).removeClass('noprint');
			$(getElementSelector(container, elementIds.ionTableLoc1)).removeClass('noprint');
			$(getElementSelector(container, elementIds.ionTableLoc2)).removeClass('noprint');
			$(getElementSelector(container, elementIds.viewOptionsDiv)).removeClass('noprint');
			$("#tempPrintDiv").siblings().removeClass("noprint");
			


			plotOptions.series.peaks.print = false; // draw the labels in the canvas
			plotAccordingToChoices(container);
			
			// move the plots back to the original location
            parent.append(container.detach());
			$("#tempPrintDiv").remove();
			
			
			/*var canvas = plot.getCanvas();
			var iWidth=3500;
			var iHeight = 3050;
			var oSaveCanvas = document.createElement("canvas");
			oSaveCanvas.width = iWidth;
			oSaveCanvas.height = iHeight;
			oSaveCanvas.style.width = iWidth+"px";
			oSaveCanvas.style.height = iHeight+"px";
			var oSaveCtx = oSaveCanvas.getContext("2d");
			oSaveCtx.drawImage(canvas, 0, 0, canvas.width, canvas.height, 0, 0, iWidth, iHeight);
			
			var dataURL = oSaveCanvas.toDataURL("image/png");
			window.location = dataURL;*/
			
			
		});
	}
	
	// -----------------------------------------------
	// SELECTED DATASETS
	// -----------------------------------------------
	function getDatasets(container) {

        var options = container.data("options");

		 // selected ions
		var selectedIonTypes = getSelectedIonTypes(container);
		
		calculateTheoreticalSeries(container, selectedIonTypes);
		
		// add the un-annotated peaks
		var data = [{data: options.peaks, color: "#bbbbbb", labelType: 'none'}];
		
		// add the annotated peaks
		var seriesMatches = getSeriesMatches(container, selectedIonTypes);
		for(var i = 0; i < seriesMatches.length; i += 1) {
			data.push(seriesMatches[i]);
		}
		
		// add any user specified extra peaks
		for(var i = 0; i < options.extraPeakSeries.length; i += 1) {
			data.push(options.extraPeakSeries[i]);
		}
		return data;
	}
	
	//-----------------------------------------------
	// SELECTED ION TYPES
	// -----------------------------------------------
	function getSelectedIonTypes(container) {

		var ions = [];
		var charges = [];
		$(getElementSelector(container, elementIds.ion_choice)).find("input:checked").each(function () {
	        var key = $(this).attr("id");
	        var tokens = key.split("_");
	        ions.push(tokens[0]);
	        charges.push(tokens[1]);
	  	});
	    
	    var selected = [];
        var ion;
        for (var i = 0; i < ions.length; i += 1) {
            selected.push(ion = Ion.get(ions[i], charges[i]));
        }
	    
	    return selected;
	}
	
	function getSelectedNtermIons(selectedIonTypes) {
		var ntermIons = [];
		
		for(var i = 0; i < selectedIonTypes.length; i += 1) {
			var sion = selectedIonTypes[i];
			if(sion.type == "a" || sion.type == "b" || sion.type == "c") 
				ntermIons.push(sion);
		}
		ntermIons.sort(function(m,n) {
			if(m.type == n.type) {
				return (m.charge - n.charge);
			}
			else {
				return m.type - n.type;
			}
		});
		return ntermIons;
	}

	function getSelectedCtermIons(selectedIonTypes) {
		var ctermIons = [];
		
		for(var i = 0; i < selectedIonTypes.length; i += 1) {
			var sion = selectedIonTypes[i];
			if(sion.type == "x" || sion.type == "y" || sion.type == "z") 
				ctermIons.push(sion);
		}
		ctermIons.sort(function(m,n) {
			if(m.type == n.type) {
				return (m.charge - n.charge);
			}
			else {
				return m.type - n.type;
			}
		});
		return ctermIons;
	}
	
	// ---------------------------------------------------------
	// CALCUALTE THEORETICAL MASSES FOR THE SELECTED ION SERIES
	// ---------------------------------------------------------
	function calculateTheoreticalSeries(container, selectedIons) {

		if(selectedIons) {
		
			var todoIonSeries = [];
			var todoIonSeriesData = [];
            var ionSeries = container.data("ionSeries");
			for(var i = 0; i < selectedIons.length; i += 1) {
				var sion = selectedIons[i];
				if(sion.type == "a") {
					if(!container.data("massTypeChanged") && ionSeries.a[sion.charge])	continue; // already calculated
					else {
						todoIonSeries.push(sion);
						ionSeries.a[sion.charge] = [];
						todoIonSeriesData.push(ionSeries.a[sion.charge]);
					}
				}
				if(sion.type == "b") {
					if(!container.data("massTypeChanged") && ionSeries.b[sion.charge])	continue; // already calculated
					else {
						todoIonSeries.push(sion);
						ionSeries.b[sion.charge] = [];
						todoIonSeriesData.push(ionSeries.b[sion.charge]);
					}
				}
				if(sion.type == "c") {
					if(!container.data("massTypeChanged") && ionSeries.c[sion.charge])	continue; // already calculated
					else {
						todoIonSeries.push(sion);
						ionSeries.c[sion.charge] = [];
						todoIonSeriesData.push(ionSeries.c[sion.charge]);
					}
				}
				if(sion.type == "x") {
					if(!container.data("massTypeChanged") && ionSeries.x[sion.charge])	continue; // already calculated
					else {
						todoIonSeries.push(sion);
						ionSeries.x[sion.charge] = [];
						todoIonSeriesData.push(ionSeries.x[sion.charge]);
					}
				}
				if(sion.type == "y") {
					if(!container.data("massTypeChanged") && ionSeries.y[sion.charge])	continue; // already calculated
					else {
						todoIonSeries.push(sion);
						ionSeries.y[sion.charge] = [];
						todoIonSeriesData.push(ionSeries.y[sion.charge]);
					}
				}
				if(sion.type == "z") {
					if(!container.data("massTypeChanged") && ionSeries.z[sion.charge])	continue; // already calculated
					else {
						todoIonSeries.push(sion);
						ionSeries.z[sion.charge] = [];
						todoIonSeriesData.push(ionSeries.z[sion.charge]);
					}
				}
			}

			if(container.data("options").sequence) {

                var sequence =  container.data("options").sequence
				var massType = container.find("input[name='"+getRadioName(container, "massTypeOpt")+"']:checked").val();
				
				for(var i = 1; i < sequence.length; i += 1) {
					
					for(var j = 0; j < todoIonSeries.length; j += 1) {
						var tion = todoIonSeries[j];
						var ionSeriesData = todoIonSeriesData[j];
						if(tion.term == "n")
							ionSeriesData.push(sion = Ion.getSeriesIon(tion, container.data("options").peptide, i, massType));
						else if(tion.term == "c")
							ionSeriesData.unshift(sion = Ion.getSeriesIon(tion, container.data("options").peptide, i, massType));
					}
				}
			}
		}
	}

	
	// -----------------------------------------------
	// MATCH THEORETICAL MASSES WITH PEAKS IN THE SCAN
	// -----------------------------------------------
	function recalculate(container) {
        return (container.data("massErrorChanged") ||
				container.data("massTypeChanged") ||
				container.data("peakAssignmentTypeChanged") ||
				container.data("selectedNeutralLossChanged"));
	}

	function getSeriesMatches(container, selectedIonTypes) {
		
		var dataSeries = [];
		
		var peakAssignmentType = container.find("input[name='"+getRadioName(container, "peakAssignOpt")+"']:checked").val();
		var peakLabelType = container.find("input[name='"+getRadioName(container, "peakLabelOpt")+"']:checked").val();

        var ionSeriesMatch = container.data("ionSeriesMatch");
        var ionSeries = container.data("ionSeries");
        var ionSeriesLabels = container.data("ionSeriesLabels");
        var options = container.data("options");
        var massError = container.data("massError");
        var peaks = getPeaks(container);

		for(var j = 0; j < selectedIonTypes.length; j += 1) {
		
			var ion = selectedIonTypes[j];
							
			
			if(ion.type == "a") {
				if(recalculate(container) || !ionSeriesMatch.a[ion.charge]) { // re-calculate only if mass error has changed OR
																		// matching peaks for this series have not been calculated
					// calculated matching peaks
					var adata = calculateMatchingPeaks(container, ionSeries.a[ion.charge], peaks, massError, peakAssignmentType);
					if(adata && adata.length > 0) {
						ionSeriesMatch.a[ion.charge] = adata[0];
						ionSeriesLabels.a[ion.charge] = adata[1];
					}
				}
				dataSeries.push({data: ionSeriesMatch.a[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.a[ion.charge]});
			}
			
			if(ion.type == "b") {
				if(recalculate(container) || !ionSeriesMatch.b[ion.charge]) { // re-calculate only if mass error has changed OR
																		// matching peaks for this series have not been calculated
					// calculated matching peaks
					var bdata = calculateMatchingPeaks(container, ionSeries.b[ion.charge], peaks, massError, peakAssignmentType);
					if(bdata && bdata.length > 0) {
						ionSeriesMatch.b[ion.charge] = bdata[0];
						ionSeriesLabels.b[ion.charge] = bdata[1];
					}
				}
				dataSeries.push({data: ionSeriesMatch.b[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.b[ion.charge]});
			}
			
			if(ion.type == "c") {
				if(recalculate(container) || !ionSeriesMatch.c[ion.charge]) { // re-calculate only if mass error has changed OR
																		// matching peaks for this series have not been calculated
					// calculated matching peaks
					var cdata = calculateMatchingPeaks(container, ionSeries.c[ion.charge], peaks, massError, peakAssignmentType);
					if(cdata && cdata.length > 0) {
						ionSeriesMatch.c[ion.charge] = cdata[0];
						ionSeriesLabels.c[ion.charge] = cdata[1];
					}
				}
				dataSeries.push({data: ionSeriesMatch.c[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.c[ion.charge]});
			}
			
			if(ion.type == "x") {
				if(recalculate(container) || !ionSeriesMatch.x[ion.charge]) { // re-calculate only if mass error has changed OR
																		// matching peaks for this series have not been calculated
					// calculated matching peaks
					var xdata = calculateMatchingPeaks(container, ionSeries.x[ion.charge], peaks, massError, peakAssignmentType);
					if(xdata && xdata.length > 0) {
						ionSeriesMatch.x[ion.charge] = xdata[0];
						ionSeriesLabels.x[ion.charge] = xdata[1];
					}
				}
				dataSeries.push({data: ionSeriesMatch.x[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.x[ion.charge]});
			}
			
			if(ion.type == "y") {
				if(recalculate(container) || !ionSeriesMatch.y[ion.charge]) { // re-calculate only if mass error has changed OR
																		// matching peaks for this series have not been calculated
					// calculated matching peaks
					var ydata = calculateMatchingPeaks(container, ionSeries.y[ion.charge], peaks, massError, peakAssignmentType);
					if(ydata && ydata.length > 0) {
						ionSeriesMatch.y[ion.charge] = ydata[0];
						ionSeriesLabels.y[ion.charge] = ydata[1];
					}
				}
				dataSeries.push({data: ionSeriesMatch.y[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.y[ion.charge]});
			}
			
			if(ion.type == "z") {
				if(recalculate(container) || !ionSeriesMatch.z[ion.charge]) { // re-calculate only if mass error has changed OR
																		// matching peaks for this series have not been calculated
					// calculated matching peaks
					var zdata = calculateMatchingPeaks(container, ionSeries.z[ion.charge], peaks, massError, peakAssignmentType);
					if(zdata && zdata.length > 0) {
						ionSeriesMatch.z[ion.charge] = zdata[0];
						ionSeriesLabels.z[ion.charge] = zdata[1];
					}
				}
				dataSeries.push({data: ionSeriesMatch.z[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.z[ion.charge]});
			}
		}
		return dataSeries;
	}

	function calculateMatchingPeaks(container, ionSeries, allPeaks, massTolerance, peakAssignmentType) {

        // console.log("calculating matching peaks");
		var peakIndex = 0;
		
		var matchData = [];
		matchData[0] = []; // peaks
		matchData[1] = []; // labels -- ions;
		
		var neutralLosses = [];
		$(getElementSelector(container, elementIds.nl_choice)).find("input:checked").each(function() {
			neutralLosses.push($(this).val());
		});
		for(var i = 0; i < ionSeries.length; i += 1) {
			
			var sion = ionSeries[i];
			
			// get match for water and or ammonia loss
			for(var n = 0; n < neutralLosses.length; n += 1) {
				getMatchForIon(sion, matchData, allPeaks, peakIndex, massTolerance, peakAssignmentType, neutralLosses[n]);
			}
			// get match for the ion
			peakIndex = getMatchForIon(sion, matchData, allPeaks, peakIndex, massTolerance, peakAssignmentType);
		}
		
		return matchData;
	}
	
	// sion -- theoretical ion
	// matchData -- array to which we will add a peak if there is a match
	// allPeaks -- array with all the scan peaks
	// peakIndex -- current index in peaks array
	// Returns the index of the matching peak, if one is found
	function getMatchForIon(sion, matchData, allPeaks, peakIndex, massTolerance, peakAssignmentType, neutralLoss) {
		
		var bestPeak = null;
		if(!neutralLoss)
			sion.match = false; // reset;
		var ionmz;
		if(!neutralLoss)
			ionmz = sion.mz;
		else {
			if(neutralLoss == 'h2o') {
				ionmz = Ion.getWaterLossMz(sion);
			}
			else if(neutralLoss = 'nh3') {
				ionmz = Ion.getAmmoniaLossMz(sion);
			}
		}
		var bestDistance;
		
		for(var j = peakIndex; j < allPeaks.length; j += 1) {
			
			var peak = allPeaks[j];
			
			// peak is before the current ion we are looking at
			if(peak[0] < ionmz - massTolerance)
				continue;
				
			// peak is beyond the current ion we are looking at
			if(peak[0] > ionmz + massTolerance) {
			
				// if we found a matching peak for the current ion, save it
				if(bestPeak) {
					//console.log("found match "+sion.label+", "+ionmz+";  peak: "+bestPeak[0]);
					matchData[0].push([bestPeak[0], bestPeak[1]]);
					if(!neutralLoss) {
						matchData[1].push(sion.label);
						sion.match = true;
					}
					else {
						if(neutralLoss == 'h2o') {
							matchData[1].push(sion.label+'o');
						}
						else if(neutralLoss = 'nh3') {
							matchData[1].push(sion.label+'*');
						}
					}
				}
				peakIndex = j;
				break;
			}
				
			// peak is within +/- massTolerance of the current ion we are looking at
			
			// if this is the first peak in the range
			if(!bestPeak) {
				//console.log("found a peak in range, "+peak.mz);
				bestPeak = peak;
				bestDistance = Math.abs(ionmz - peak[0]);
				continue;
			}
			
			// if peak assignment method is Most Intense
			if(peakAssignmentType == "intense") {
				if(peak[1] > bestPeak[1]) {
					bestPeak = peak;
					continue;
				}
			}
			
			// if peak assignment method is Closest Peak
			if(peakAssignmentType == "close") {
				var dist = Math.abs(ionmz - peak[0]);
				if(!bestDistance || dist < bestDistance) {
					bestPeak = peak;
					bestDistance = dist;
				}
			}
		}
		
		return peakIndex;
	}
	

    function getPeaks(container)
    {
        var options = container.data("options");

        if($(getElementSelector(container, elementIds.removeNoise)).is(":checked"))
        {
            if(options.sparsePeaks == null) {
                calculateSparsePeaks(container);
            }
            return options.sparsePeaks;
        }
        else
        {
            return options.peaks;
        }
    }

    function calculateSparsePeaks(container) {

        console.log("calculating sparse peaks");

        var peaks = container.data("options").peaks;
        var sparsePeaks = [];

        for(var i = 0; i < peaks.length; i += 1) {

			var peak = peaks[i];

            var intensity = peak[1];
            var mz = peak[0];
            var minMz = mz;
            var maxMz = mz;
            var j = i-1;
            var totalIntensity = intensity;
            var peakCount = 1;
            // sum up the intensities in the +/- 50Da window of this peak
            while((minMz >= mz - 50.0) && j >= 0)
            {
                totalIntensity += peaks[j][1];
                minMz = peaks[j][0];
                j -= 1;
                peakCount += 1;
            }
            j = i+1;
            while(maxMz <= mz + 50.0 && j < peaks.length)
            {
                totalIntensity += peaks[j][1];
                maxMz = peaks[j][0];
                j += 1;
                peakCount += 1;
            }

            var mean = totalIntensity / peakCount;

            // calculate the standard deviation
            var sdev = 0;
            j = i - 1;
            while((minMz >= mz - 50.0) && j >= 0)
            {
                sdev += Math.pow((peaks[j][1] - mean), 2);
                minMz = peaks[j][0];
                j -= 1;
            }
            j = i+1;
            while(maxMz <= mz + 50.0 && j < peaks.length)
            {
                sdev += Math.pow((peaks[j][1] - mean), 2);
                maxMz = peaks[j][0];
                j += 1;
            }
            sdev = Math.sqrt(sdev / peakCount);

            if(intensity >= mean + 2 * sdev)
            {
                sparsePeaks.push(peak);
            }
            //console.log(intensity+"  "+mean+"  "+sdev);
		}
        console.log("Sparse Peak count: "+sparsePeaks.length);
        console.log("All Peaks count: "+peaks.length);
        container.data("options").sparsePeaks = sparsePeaks;
    }

	// -----------------------------------------------
	// INITIALIZE THE CONTAINER
	// -----------------------------------------------
    function createContainer(div) {

        div.append('<div id="'+elementIds.lorikeet_content+"_"+index+'"></div>');
        var container = $("#"+ div.attr('id')+" > #"+elementIds.lorikeet_content+"_"+index);
        container.addClass("lorikeet");
        return container;
    }

	function initContainer(container) {

        var options = container.data("options");

        var rowspan = 2;

		var parentTable = '<table cellpadding="0" cellspacing="5"> ';
		parentTable += '<tbody> ';
		parentTable += '<tr> ';

		// Header
		parentTable += '<td colspan="3" class="bar"> ';
		parentTable += '</div> ';
		parentTable += '</td> ';
		parentTable += '</tr> ';

		// options table
		parentTable += '<tr> ';
		parentTable += '<td rowspan="'+rowspan+'" valign="top" id="'+getElementId(container, elementIds.optionsTable)+'"> ';
		parentTable += '</td> ';

		// placeholder for sequence, m/z, scan number etc
		parentTable += '<td style="background-color: white; padding:5px; border:1px dotted #cccccc;" valign="bottom" align="center"> ';
		parentTable += '<div id="'+getElementId(container, elementIds.seqinfo)+'" style="width:100%;"></div> ';
		// placeholder for file name, scan number and charge
		parentTable += '<div id="'+getElementId(container, elementIds.fileinfo)+'" style="width:100%;"></div> ';
		parentTable += '</td> ';


		// placeholder for the ion table
		parentTable += '<td rowspan="'+rowspan+'" valign="top" id="'+getElementId(container, elementIds.ionTableLoc1)+'" > ';
		parentTable += '<div id="'+getElementId(container, elementIds.ionTableDiv)+'">';
		parentTable += '<span id="'+getElementId(container, elementIds.moveIonTable)+'" class="font_small link">[Click]</span> <span class="font_small">to move table</span>';
		// placeholder for file name, scan number, modifications etc.
		parentTable += '<div id="'+getElementId(container, elementIds.modInfo)+'" style="margin-top:5px;"></div> ';
		parentTable += '</div> ';
		parentTable += '</td> ';
		parentTable += '</tr> ';


		// placeholders for the ms/ms plot
		parentTable += '<tr> ';
		parentTable += '<td style="background-color: white; padding:5px; border:1px dotted #cccccc;" valign="middle" align="center"> ';
		parentTable += '<div id="'+getElementId(container, elementIds.msmsplot)+'" align="bottom" style="width:'+options.width+'px;height:'+options.height+'px;"></div> ';

		// placeholder for viewing options (zoom, plot size etc.)
		parentTable += '<div id="'+getElementId(container, elementIds.viewOptionsDiv)+'" align="top" style="margin-top:15px;"></div> ';

		// placeholder for ms1 plot (if data is available)
		if(options.ms1peaks && options.ms1peaks.length > 0) {
			parentTable += '<div id="'+getElementId(container, elementIds.msPlot)+'" style="width:'+options.width+'px;height:100px;"></div> ';
		}
		parentTable += '</td> ';
		parentTable += '</tr> ';


		// Footer & placeholder for moving ion table
		parentTable += '<tr> ';
		parentTable += '<td colspan="3" class="bar noprint" valign="top" align="center" id="'+getElementId(container, elementIds.ionTableLoc2)+'" > ';
		parentTable += '<div align="center" style="width:100%;font-size:10pt;"> ';
		parentTable += '</div> ';
		parentTable += '</td> ';
		parentTable += '</tr> ';

		parentTable += '</tbody> ';
		parentTable += '</table> ';

		container.append(parentTable);
		
		return container;
	}
	
	
	//---------------------------------------------------------
	// ION TABLE
	//---------------------------------------------------------
	function makeIonTable(container) {

        var options = container.data("options");

	 	// selected ions
		var selectedIonTypes = getSelectedIonTypes(container);
		var ntermIons = getSelectedNtermIons(selectedIonTypes);
		var ctermIons = getSelectedCtermIons(selectedIonTypes);
		
		var myTable = '' ;
		myTable += '<table id="'+getElementId(container, elementIds.ionTable)+'" cellpadding="2" class="font_small '+elementIds.ionTable+'">' ;
		myTable +=  "<thead>" ;
		myTable +=   "<tr>";

		// nterm ions
		for(var i = 0; i < ntermIons.length; i += 1) {
			myTable +=    "<th>" +ntermIons[i].label+  "</th>";   
		}
		myTable +=    "<th>" +"#"+  "</th>"; 
		myTable +=    "<th>" +"Seq"+  "</th>"; 
		myTable +=    "<th>" +"#"+  "</th>"; 
		// cterm ions
		for(var i = 0; i < ctermIons.length; i += 1) {
			myTable +=    "<th>" +ctermIons[i].label+  "</th>"; 
		}
		myTable +=   "</tr>" ;
		myTable +=  "</thead>" ;
		
		myTable +=  "<tbody>" ;

        var ionSeries = container.data("ionSeries");

		for(var i = 0; i < options.sequence.length; i += 1) {
            var aaChar = options.sequence.charAt(i);
			myTable +=   "<tr>";

			// nterm ions
			for(var n = 0; n < ntermIons.length; n += 1) {
				if(i < options.sequence.length - 1) {
					var seriesData = getCalculatedSeries(ionSeries, ntermIons[n]);
					var cls = "";
					var style = "";
					if(seriesData[i].match) {
						cls="matchIon";
						style="style='background-color:"+Ion.getSeriesColor(ntermIons[n])+";'";
					}
					myTable +=    "<td class='"+cls+"' "+style+" >" +round(seriesData[i].mz)+  "</td>";  
				}
				else {
					myTable +=    "<td>" +"&nbsp;"+  "</td>"; 
				} 
			}
			
			myTable += "<td class='numCell'>"+(i+1)+"</td>";
			if(options.peptide.varMods()[i+1])
				myTable += "<td class='seq modified'>"+aaChar+"</td>";
			else
				myTable += "<td class='seq'>"+aaChar+"</td>";
			myTable += "<td class='numCell'>"+(options.sequence.length - i)+"</td>";
			
			// cterm ions
			for(var c = 0; c < ctermIons.length; c += 1) {
				if(i > 0) {
					var seriesData = getCalculatedSeries(ionSeries, ctermIons[c]);
					var idx = options.sequence.length - i - 1;
					var cls = "";
					var style = "";
					if(seriesData[idx].match) {
						cls="matchIon";
						style="style='background-color:"+Ion.getSeriesColor(ctermIons[c])+";'";
					}
					myTable +=    "<td class='"+cls+"' "+style+" >" +round(seriesData[idx].mz)+  "</td>";  
				}
				else {
					myTable +=    "<td>" +"&nbsp;"+  "</td>"; 
				} 
			}
			
		}
		myTable +=   "</tr>" ;
		
		myTable += "</tbody>";
		myTable += "</table>";
		
		// alert(myTable);
		$(getElementSelector(container, elementIds.ionTable)).remove();
		$(getElementSelector(container, elementIds.ionTableDiv)).prepend(myTable); // add as the first child
	}
	
	function getCalculatedSeries(ionSeries, ion) {
        if(ion.type == "a")
			return ionSeries.a[ion.charge];
		if(ion.type == "b")
			return ionSeries.b[ion.charge];
		if(ion.type == "c")
			return ionSeries.c[ion.charge];
		if(ion.type == "x")
			return ionSeries.x[ion.charge];
		if(ion.type == "y")
			return ionSeries.y[ion.charge];
		if(ion.type == "z")
			return ionSeries.z[ion.charge];
	}
	
	function makeIonTableMovable(container) {

		$(getElementSelector(container, elementIds.moveIonTable)).hover(
			function(){
		         $(this).css({cursor:'pointer'}); //mouseover
		    }
		);

		$(getElementSelector(container, elementIds.moveIonTable)).click(function() {
			var ionTableDiv = $(getElementSelector(container, elementIds.ionTableDiv));
			if(ionTableDiv.is(".moved")) {
				ionTableDiv.removeClass("moved");
				ionTableDiv.detach();
				$(getElementSelector(container, elementIds.ionTableLoc1)).append(ionTableDiv);
			}
			else {
				ionTableDiv.addClass("moved");
				ionTableDiv.detach();
				$(getElementSelector(container, elementIds.ionTableLoc2)).append(ionTableDiv);
			}
		});
	}

	//---------------------------------------------------------
	// SEQUENCE INFO
	//---------------------------------------------------------
	function showSequenceInfo (container) {

        var options = container.data("options");

		var specinfo = '';
		if(options.sequence) {
			
			specinfo += '<div>';
			specinfo += '<span style="font-weight:bold; color:#8B0000;">'+getModifiedSequence(options)+'</span>';
			
			var massType = container.find("input[name='"+getRadioName(container, "massTypeOpt")+"']:checked").val();
			
			var neutralMass = 0;
			
			if(massType == "mono")
				neutralMass = options.peptide.getNeutralMassMono();
			else if(massType == "avg")
				neutralMass = options.peptide.getNeutralMassAvg();
				
			
			var mz;
			if(options.charge) {
				mz = Ion.getMz(neutralMass, options.charge);
			}
			
			var mass = neutralMass + Ion.MASS_PROTON;
			specinfo += ', MH+ '+mass.toFixed(4);
			if(mz) 
				specinfo += ', m/z '+mz.toFixed(4);
			specinfo += '</div>';
			
		}
		
		// first clear the div if it has anything
		$(getElementSelector(container, elementIds.seqinfo)).empty();
		$(getElementSelector(container, elementIds.seqinfo)).append(specinfo);
	}
	
	function getModifiedSequence(options) {
		
		var modSeq = '';
		for(var i = 0; i < options.sequence.length; i += 1) {
			
			if(options.peptide.varMods()[i+1])
				modSeq += '<span style="background-color:yellow;padding:1px;border:1px dotted #CFCFCF;">'+options.sequence.charAt(i)+"</span>";
			else
				modSeq += options.sequence.charAt(i);
		}
		
		return modSeq;
	}
	
	//---------------------------------------------------------
	// FILE INFO -- filename, scan number, precursor m/z and charge
	//---------------------------------------------------------
	function showFileInfo (container) {

        var options = container.data("options");

		var fileinfo = '';
			
		if(options.fileName || options.scanNum) {
			fileinfo += '<div style="margin-top:5px;" class="font_small">';
			if(options.fileName) {
				fileinfo += 'File: '+options.fileName;
			}
			if(options.scanNum) {
				fileinfo += ', Scan: '+options.scanNum;
			}	
			if(options.precursorMz) {
				fileinfo += ', Precursor m/z: '+options.precursorMz;
			}
			if(options.charge) {
				fileinfo += ', Charge: '+options.charge;
			}
			fileinfo += '</div>';
		}
		
		$(getElementSelector(container, elementIds.fileinfo)).append(fileinfo);
	}
	
	//---------------------------------------------------------
	// MODIFICATION INFO
	//---------------------------------------------------------
	function showModInfo (container) {

        var options = container.data("options");

		var modInfo = '';
			
		modInfo += '<div>';
		if(options.ntermMod || options.ntermMod > 0) {
			modInfo += 'Add to N-term: <b>'+options.ntermMod+'</b>';
		}
		if(options.ctermMod || options.ctermMod > 0) {
			modInfo += 'Add to C-term: <b>'+options.ctermMod+'</b>';
		}
		modInfo += '</div>';
		
		if(options.staticMods && options.staticMods.length > 0) {
			modInfo += '<div style="margin-top:5px;">';
			modInfo += 'Static Modifications: ';
			for(var i = 0; i < options.staticMods.length; i += 1) {
				var mod = options.staticMods[i];
				//if(i > 0) modInfo += ', ';
				modInfo += "<div><b>"+mod.aa.code+": "+mod.modMass+"</b></div>";
			}
			modInfo += '</div>';
		}
		
		if(options.variableMods && options.variableMods.length > 0) {
			
			var modChars = [];
			var uniqvarmods = [];
			for(var i = 0; i < options.variableMods.length; i += 1) {
				var mod = options.variableMods[i];
				if(modChars[mod.aa.code])
					continue;
				modChars[mod.aa.code] = 1;
				uniqvarmods.push(mod);
			}  
			
			modInfo += '<div style="margin-top:5px;">';
			modInfo += 'Variable Modifications: ';
			for(var i = 0; i < uniqvarmods.length; i += 1) {
				var mod = uniqvarmods[i];
				//if(i > 0) modInfo += ', ';
				modInfo += "<div><b>"+mod.aa.code+": "+mod.modMass+"</b></div>";
			}
			modInfo += '</div>';
		}
		
		$(getElementSelector(container, elementIds.modInfo)).append(modInfo);
	}
	
	//---------------------------------------------------------
	// VIEWING OPTIONS TABLE
	//---------------------------------------------------------
	function makeViewingOptions(container) {

        var options = container.data("options");

		var myContent = '';
		
		// reset zoom option
		myContent += '<nobr> ';
		myContent += '<span style="width:100%; font-size:8pt; margin-top:5px; color:sienna;">Click and drag in the plot to zoom</span> ';
		myContent += 'X:<input id="'+getElementId(container, elementIds.zoom_x)+'" type="checkbox" value="X" checked="checked"/> ';
		myContent += '&nbsp;Y:<input id="'+getElementId(container, elementIds.zoom_y)+'" type="checkbox" value="Y" /> ';
		myContent += '&nbsp;<input id="'+getElementId(container, elementIds.resetZoom)+'" type="button" value="Zoom Out" /> ';
		myContent += '&nbsp;<input id="'+getElementId(container, elementIds.printLink)+'" type="button" value="Print" /> ';
		myContent += '</nobr> ';
		
		myContent += '&nbsp;&nbsp;';
		
		// tooltip option
		myContent += '<nobr> ';
		myContent += '<input id="'+getElementId(container, elementIds.enableTooltip)+'" type="checkbox">Enable tooltip ';
		myContent += '</nobr> ';
		
		myContent += '<br>';
		
		$(getElementSelector(container, elementIds.viewOptionsDiv)).append(myContent);
		if(!options.showViewingOptions) {
            $(getElementSelector(container, elementIds.viewOptionsDiv)).hide();
        }
	}
	
	
	//---------------------------------------------------------
	// OPTIONS TABLE
	//---------------------------------------------------------
	function makeOptionsTable(container, charge) {

        var options = container.data("options");

		var myTable = '';
		myTable += '<table cellpadding="2" cellspacing="2"> ';
		myTable += '<tbody> ';
		
		// Ions
		myTable += '<tr><td class="optionCell"> ';
		
		myTable += '<b>Ions:</b> ';
		myTable += '<div id="'+getElementId(container, elementIds.ion_choice)+'" style="margin-bottom: 10px"> ';
		myTable += '<!-- a ions --> ';
		myTable += '<nobr> ';
		myTable += '<span style="font-weight: bold;">a</span> ';
		myTable += '<input type="checkbox" value="1" id="a_1"/>1<sup>+</sup> ';
		myTable += '<input type="checkbox" value="2" id="a_2"/>2<sup>+</sup> ';
		myTable += '<input type="checkbox" value="3" id="a_3"/>3<sup>+</sup> ';
		myTable += '</nobr> ';
		myTable += '<br/> ';
		myTable += '<!-- b ions --> ';
		myTable += '<nobr> ';
		myTable += '<span style="font-weight: bold;">b</span> ';
		myTable += '<input type="checkbox" value="1" id="b_1" checked="checked"/>1<sup>+</sup> ';
		if(charge >= 2)
			myTable += '<input type="checkbox" value="2" id="b_2" checked="checked"/>2<sup>+</sup> ';
		else
			myTable += '<input type="checkbox" value="2" id="b_2"/>2<sup>+</sup> ';
		if(charge >= 3)
			myTable += '<input type="checkbox" value="3" id="b_3" checked="checked"/>3<sup>+</sup> ';
		else
			myTable += '<input type="checkbox" value="3" id="b_3"/>3<sup>+</sup> ';
		myTable += '</nobr> ';
		myTable += '<br/> ';
		myTable += '<!-- c ions --> ';
		myTable += '<nobr> ';
		myTable += '<span style="font-weight: bold;">c</span> ';
		myTable += '<input type="checkbox" value="1" id="c_1"/>1<sup>+</sup> ';
		myTable += '<input type="checkbox" value="2" id="c_2"/>2<sup>+</sup> ';
		myTable += '<input type="checkbox" value="3" id="c_3"/>3<sup>+</sup> ';
		myTable += '</nobr> ';
		myTable += '<br/> ';
		myTable += '<!-- x ions --> ';
		myTable += '<nobr> ';
		myTable += '<span style="font-weight: bold;">x</span> ';
		myTable += '<input type="checkbox" value="1" id="x_1"/>1<sup>+</sup> ';
		myTable += '<input type="checkbox" value="2" id="x_2"/>2<sup>+</sup> ';
		myTable += '<input type="checkbox" value="3" id="x_3"/>3<sup>+</sup> ';
		myTable += '</nobr> ';
		myTable += '<br/> ';
		myTable += '<!-- y ions --> ';
		myTable += '<nobr> ';
		myTable += '<span style="font-weight: bold;">y</span> ';
		myTable += '<input type="checkbox" value="1" id="y_1" checked="checked"/>1<sup>+</sup> ';
		if(charge >= 2)
			myTable += '<input type="checkbox" value="2" id="y_2" checked="checked"/>2<sup>+</sup> ';
		else 
			myTable += '<input type="checkbox" value="2" id="y_2"/>2<sup>+</sup> ';
		if(charge >= 3)
			myTable += '<input type="checkbox" value="3" id="y_3" checked="checked"/>3<sup>+</sup> ';
		else
			myTable += '<input type="checkbox" value="3" id="y_3"/>3<sup>+</sup> ';
		myTable += '</nobr> ';
		myTable += '<br/> ';
		myTable += '<!-- z ions --> ';
		myTable += '<nobr> ';
		myTable += '<span style="font-weight: bold;">z</span> ';
		myTable += '<input type="checkbox" value="1" id="z_1"/>1<sup>+</sup> ';
		myTable += '<input type="checkbox" value="2" id="z_2"/>2<sup>+</sup> ';
		myTable += '<input type="checkbox" value="3" id="z_3"/>3<sup>+</sup> ';
		myTable += '</nobr> ';
		myTable += '<br/> ';
		myTable += '<span id="'+getElementId(container, elementIds.deselectIonsLink)+'" style="font-size:8pt;text-decoration: underline; color:sienna;cursor:pointer;">[Deselect All]</span> ';
		myTable += '</div> ';
		
		myTable += '<span style="font-weight: bold;">Neutral Loss:</span> ';
		myTable += '<div id="'+getElementId(container, elementIds.nl_choice)+'"> ';
		myTable += '<nobr> <input type="checkbox" value="h2o" id="h2o"/> ';
		myTable += ' H<sub>2</sub>O (<span style="font-weight: bold;">o</span>)';
		myTable += '</nobr> ';
		myTable += '<br> ';
		myTable += '<nobr> <input type="checkbox" value="nh3" id="nh3"/> ';
		myTable += ' NH<sub>3</sub> (<span style="font-weight: bold;">*</span>)';
		myTable += '</nobr> ';
		myTable += '</div> ';
		
		myTable += '</td> </tr> ';
		
		// mass type, mass tolerance etc.
		myTable += '<tr><td class="optionCell"> ';
		myTable += '<div> Mass Type:<br/> ';
		myTable += '<nobr> ';
		myTable += '<input type="radio" name="'+getRadioName(container, "massTypeOpt")+'" value="mono" checked="checked"/><span style="font-weight: bold;">Mono</span> ';
		myTable += '<input type="radio" name="'+getRadioName(container, "massTypeOpt")+'" value="avg"/><span style="font-weight: bold;">Avg</span> ';
		myTable += '</nobr> ';
		myTable += '</div> ';
		myTable += '<div style="margin-top:10px;"> ';
		myTable += '<nobr>Mass Tol: <input id="'+getElementId(container, elementIds.massError)+'" type="text" value="'+options.massError+'" size="4"/></nobr> ';
		myTable += '</div> ';
		myTable += '<div style="margin-top:10px;" align="center"> ';
		myTable += '<input id="'+getElementId(container, elementIds.update)+'" type="button" value="Update"/> ';
		myTable += '</div> ';
		myTable += '</td> </tr> ';
		
		// peak assignment method
		myTable += '<tr><td class="optionCell"> ';
		myTable+= '<div> Peak Assignment:<br/> ';
		myTable+= '<input type="radio" name="'+getRadioName(container, "peakAssignOpt")+'" value="intense" checked="checked"/><span style="font-weight: bold;">Most Intense</span><br/> ';
		myTable+= '<input type="radio" name="'+getRadioName(container, "peakAssignOpt")+'" value="close"/><span style="font-weight: bold;">Nearest Match</span><br/> ';
        myTable+= '<input type="checkbox" value="true" checked="checked" id="'+getElementId(container, elementIds.removeNoise)+'"/><span style="font-weight:bold;">Remove Noise</span>';
		myTable+= '</div> ';
		myTable += '</td> </tr> ';
		
		// peak labels
		myTable += '<tr><td class="optionCell"> ';
		myTable+= '<div> Peak Labels:<br/> ';
		myTable+= '<input type="radio" name="'+getRadioName(container, "peakLabelOpt")+'" value="ion" checked="checked"/><span style="font-weight: bold;">Ion</span>';
		myTable+= '<input type="radio" name="'+getRadioName(container, "peakLabelOpt")+'" value="mz"/><span style="font-weight: bold;">m/z</span><br/>';
		myTable+= '<input type="radio" name="'+getRadioName(container, "peakLabelOpt")+'" value="none"/><span style="font-weight: bold;">None</span> ';
		myTable+= '</div> ';
		myTable += '</td> </tr> ';
		
		// sliders to change plot size
		myTable += '<tr><td class="optionCell"> ';
		myTable += '<div>Width: <span id="'+getElementId(container, elementIds.slider_width_val)+'">'+options.width+'</span></div> ';
		myTable += '<div id="'+getElementId(container, elementIds.slider_width)+'" style="margin-bottom:15px;"></div> ';
		myTable += '<div>Height: <span id="' + getElementId(container, elementIds.slider_height_val) + '">' + options.height + '</span></div> ';
		myTable += '<div id="'+getElementId(container, elementIds.slider_height)+'"></div> ';
		myTable += '</td> </tr> ';

		myTable += '</tbody>';
		myTable += '</table>';

		$(getElementSelector(container, elementIds.optionsTable)).append(myTable);
        if(!options.showIonTable) {
            $(getElementSelector(container, elementIds.optionsTable)).hide();
        }
	}

	
})(jQuery);