// $LastChangedDate$
// $LastChangedBy$
// $LastChangedRevision$

(function($) {

    // plugin name - specview
    $.fn.specview = function (opts) {

        var defaults = {
                sequence: null,
                scanNum: null,
                fileName: null,
                charge: null,
                fragmentMassType: 'mono',
                precursorMassType: 'mono',
                peakDetect: true,
                calculatedMz: null,
                precursorMz: null,
                selWinLow: null,
                selWinHigh: null,
                precursorIntensity: null,
                staticMods: [],
                variableMods: [],
                ntermMod: 0, // additional mass to be added to the n-term
                ctermMod: 0, // additional mass to be added to the c-term
                maxNeutralLossCount: 1,
                peaks: [],
                showA:[],
                showB:[],
                showC:[],
                showX:[],
                showY:[],
                showZ:[],
                sparsePeaks: null,
                ms1peaks: null,
                ms1scanLabel: null,
                ms2peaks2Label: '',
                precursorPeaks: null,
                precursorPeakClickFn: null,
                zoomMs1: false,
                width: 700, 	// width of the ms/ms plot
                height: 450, 	// height of the ms/ms plot
                massError: 0.5, // mass tolerance (in th) for labeling peaks
                massErrorUnit: massErrorTypeTh, // 'Th' or 'ppm'
                extraPeakSeries:[],
                showIonTable: true,
                showViewingOptions: true,
                showOptionsTable: true,
                peakLabelOpt: 'mz',
                showSequenceInfo: true,
                labelImmoniumIons: true,
                labelPrecursorPeak: true,
                labelReporters: false,
	        tooltipZIndex: null,
                showMassErrorPlot: false,
                massErrorPlotDefaultUnit: null,
	        userReporterIons: null,
                // Use these options to set the x-axis range (m/z) that will be displayed when the MS/MS plot is initialized or is fully zoomed out.
                // Default range is the min and max peak m/z.
                minDisplayMz: null,
	        maxDisplayMz: null
        };

	var options = $.extend(true, {}, defaults, opts); // this is a deep copy

        return this.each(function() {

            index = index + 1;
            init($(this), options);

        });
    };

    var index = 0;
    var massErrorTypeTh = 'Th';		// Th are units of m/z
    var massErrorTypePpm = 'ppm';

    var elementIds = {
            massError: "massError",
            massErrorUnit: "massErrorUnit",
            msPlot: "msPlot",
            massErrorPlot: "massErrorPlot",
            massErrorPlot_unit: "massErrorPlot_unit",
            massErrorPlot_option: "massErrorPlot_option",
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
            peakDetect: "peakDetect",
	    labelPrecursor: "labelPrecursor",
            immoniumIons: "immoniumIons",
	    reporterIons: "reporterIons",
	    showIonTable: "showIonTable",
	    anticInfo: "anticInfo",
	    userReporterIons: "userReporterIons"
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

        // trim any 0 intensity peaks from the end of the peaks array
        trimPeaksArray(options);

        // read the static modifications
        var parsedStaticMods = [];
        for(var i = 0; i < options.staticMods.length; i += 1) {
            var mod = options.staticMods[i];
            parsedStaticMods[i] = new Modification(AminoAcid.get(mod.aminoAcid), mod.modMass, mod.losses);
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
                                    AminoAcid.get(mod.aminoAcid),
                                    mod.losses
                                );
        }
        options.variableMods = parsedVarMods;

        var peptide = new Peptide(options.sequence, options.staticMods, options.variableMods,
                                  options.ntermMod, options.ctermMod, options.maxNeutralLossCount);
        options.peptide = peptide;

        // Calculate a theoretical m/z from the given sequence and charge
        if(options.sequence && options.charge) {
            var mass = options.peptide.getNeutralMassMono();
            options.calculatedMz = Ion.getMz(mass, options.charge);
        }

        if(!options.minDisplayMz)
        {
            options.minDisplayMz = options.peaks[0][0];
        }
        if(!options.maxDisplayMz)
        {
            options.maxDisplayMz = options.peaks[options.peaks.length - 1][0];
        }

        if(options.massErrorUnit !== massErrorTypeTh && options.massErrorUnit !== massErrorTypePpm)
        {
            console.warn("Invalid mass error unit given: " + options.massErrorUnit + ". Setting to " + massErrorTypeTh + ".");
            options.massErrorUnit = massErrorTypeTh;
        }
        if(!options.massErrorPlotDefaultUnit)
        {
            options.massErrorPlotDefaultUnit = options.massErrorUnit;
        }

        var container = createContainer(parent_container);
        // alert(container.attr('id')+" parent "+container.parent().attr('id'));
        storeContainerData(container, options);
        initContainer(container);

        var defaultSelectedIons = getDefaultSelectedIons(options);
        makeOptionsTable(container,[1,2,3], defaultSelectedIons);

        makeViewingOptions(container, options);

        if(options.showSequenceInfo) {
            showSequenceInfo(container, options);
            showFileInfo(container, options);
            showModInfo(container, options);
	    showUserReporterIonsInfo(container, options);
        }

        // calculate precursor ions
        calculatePrecursorPeak(container);

        // calculate immonium ions
        calculateImmoniumIons(container);

        // Calculate reporter ion labels, if required
        calculateReporterIons(container);

	// calculate user-supplied reporter ions
        if(options.userReporterIons) {
            calculateUserReporters(options, container);
        }

        createPlot(container, getDatasets(container)); // Initial MS/MS Plot

        if(options.ms1peaks && options.ms1peaks.length > 0) {

            var precursorMz = options.precursorMz;

            if(precursorMz)
            {
                // Find an actual peak closest to the precursor
                var diff = 5.0;
                var x, y;

                for(var i = 0; i < options.ms1peaks.length; i += 1) {
                    var pk = options.ms1peaks[i];
                    var d = Math.abs(pk[0] - precursorMz);
                    if(!diff || d < diff) {
                        x = pk[0];
                        y = pk[1];
                        diff = d;
                    }
                }
                if(diff <= 0.5) {
                    options.precursorIntensity = y;
                    if(!options.precursorPeaks)
                    {
                        options.precursorPeaks = [];
                    }
                    options.precursorPeaks.push([x,y]);
                }

                // Determine a zoom range
                if(options.zoomMs1)
                {
                    var maxIntensityInRange = 0;

                    for(var j = 0; j < options.ms1peaks.length; j += 1) {
                        var pk = options.ms1peaks[j];

                        if(pk[0] < options.precursorMz - 5.0)
                            continue;
                        if(pk[0] > options.precursorMz + 5.0)
                            break;
                        if(pk[1] > maxIntensityInRange)
                            maxIntensityInRange = pk[1];
                    }

                    options.maxIntensityInMs1ZoomRange = maxIntensityInRange;

                    // Set the zoom range
                    var ms1zoomRange = container.data("ms1zoomRange");
                    ms1zoomRange = {xaxis: {}, yaxis: {}};
                    ms1zoomRange.xaxis.from = options.precursorMz - 5.0;
                    ms1zoomRange.xaxis.to = options.precursorMz + 5.0;

                    ms1zoomRange.yaxis.from = 0.0;
                    ms1zoomRange.yaxis.to = options.maxIntensityInMs1ZoomRange;
                    container.data("ms1zoomRange", ms1zoomRange);
                }
            }

            createMs1Plot(container);
            setupMs1PlotInteractions(container);
        }

        setupInteractions(container, options);

        if(options.showIonTable) {
            makeIonTable(container);
        }
    }

    function getDefaultSelectedIons(options)
    {
        var userDefined = false;
        var defaultSelected = {};
        if(options.showA.length > 0)
        {
            userDefined = true;
            defaultSelected['a'] = options.showA;
        }
        if(options.showB.length > 0)
        {
            userDefined = true;
            defaultSelected['b'] = options.showB;
        }
        if(options.showC.length > 0)
        {
            userDefined = true;
            defaultSelected['c'] = options.showC;
	}
        if(options.showX.length > 0)
        {
            userDefined = true;
            defaultSelected['x'] = options.showX;
	}
        if(options.showY.length > 0)
        {
            userDefined = true;
            defaultSelected['y'] = options.showY;
        }
        if(options.showZ.length > 0)
        {
            userDefined = true;
            defaultSelected['z'] = options.showZ;
        }
        if(!userDefined)
        {
            var selected = [1];
            if(options.charge)
            {
                for (var i = 2; i < options.charge; i += 1)
                {
                    selected.push(1);
                }
            }
            defaultSelected['b'] = selected;
            defaultSelected['y'] = selected;
        }
        return defaultSelected;
    }
    // trim any 0 intensity peaks from the end of the ms/ms peaks array
    function trimPeaksArray(options)
    {
        var peaksLength = options.peaks.length;
        var lastNonZeroIntensityPeakIndex = peaksLength - 1;
        for(var i = peaksLength - 1; i >= 0; i--)
        {
            if(options.peaks[i][1] != 0.0)
            {
                lastNonZeroIntensityPeakIndex = i;
                break;
            }
        }
        if(lastNonZeroIntensityPeakIndex < peaksLength - 1)
        {
            options.peaks.splice(lastNonZeroIntensityPeakIndex+1, peaksLength - lastNonZeroIntensityPeakIndex);
        }
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
        container.data("ionSeries", {a: [], b: [], c: [], x: [], y: [], z: []}); // theoretical ion series
        container.data("ionSeriesLabels", {a: [], b: [], c: [], x: [], y: [], z: []});
        container.data("ionSeriesMatch", {a: [], b: [], c: [], x: [], y: [], z: []});
        container.data("ionSeriesAntic", {a: [], b: [], c: [], x: [], y: [], z: []}); // total annotated ion current
        container.data("anticPrecursor", 0); // total annotated ion current: precursor peak
        container.data("anticImmonium", 0); // total annotated ion current: immonium ions
        container.data("anticReporter", 0); // total annotated ion current: reporter ions
        container.data("totalIonCurrent", getTotalIonCurrent(options.peaks)); // total ion current

        var maxInt = getMaxInt(options.peaks);
        var maxIntUp = maxInt;      //get max peaks value in up peak list
        var maxIntDown = 0;         //for 2nd "butterfly" spectrum
        if(options.peaks2){
            maxIntDown = getMaxInt(options.peaks2);   //get max peaks value in down peak list
        }

	var __xrange = getPlotXRange(options);


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
            xaxis: { tickLength: 3, tickColor: "#000",
                     min: __xrange.xmin,
                     max: __xrange.xmax},
	}

        var yaxis = { tickLength: 0, tickColor: "#000",
		  max: maxInt*1.1,
		  ticks: [0, maxInt*0.1, maxInt*0.2, maxInt*0.3, maxInt*0.4, maxInt*0.5,
			  maxInt*0.6, maxInt*0.7, maxInt*0.8, maxInt*0.9, maxInt],
		  tickFormatter: function(val, axis) {return Math.round((val * 100)/maxInt)+"%";}
		}

        var yaxes = [
            { tickLength: 0, tickColor: "#000",
              max: maxIntUp*1.1,
              min: maxIntUp*-1.1,
              ticks: [0, maxIntUp*0.1, maxIntUp*0.2, maxIntUp*0.3, maxIntUp*0.4, maxIntUp*0.5,
                      maxIntUp*0.6, maxIntUp*0.7, maxIntUp*0.8, maxIntUp*0.9, maxIntUp],
              tickFormatter: function(val, axis) {return Math.round((val * 100)/maxIntUp)+"%";}
            },
            { tickLength: 0, tickColor: "#000",
              min: maxIntDown*-1.1,
              max: maxIntDown*1.1,
              position:"bottom",
              ticks: [0, maxIntDown*0.1, maxIntDown*0.2, maxIntDown*0.3, maxIntDown*0.4, maxIntDown*0.5,
                      maxIntDown*0.6, maxIntDown*0.7, maxIntDown*0.8, maxIntDown*0.9, maxIntDown],
              tickFormatter: function(val, axis) {return Math.round((val * 100)/maxIntDown)+"%";},
              transform: function (v) { return -v; },
            }
        ]

        //added double yaxis for butterfly comparing:
        if (maxIntDown > 0) {
            plotOptions['yaxes'] =  yaxes;
            plotOptions.grid.markings = [{yaxis:{from:0, to:0}, color:"#555555", lineWidth:0.5}];
            plotOptions['legend'] = {position: "se"};
	}
        else {
            plotOptions['yaxis'] =  yaxis;
	}

        container.data("plotOptions", plotOptions);
        container.data("maxInt", maxInt);
    }

    // total ion current = sum all peak intensities
    function getTotalIonCurrent(peaks) {
        var totCurrent = 0;
        for(var peak of peaks)
            totCurrent += peak[1];
        //console.log("Ion Current: "+totCurrent);
        return totCurrent;
    }

    // extract max from peak list
    function getMaxInt(peaks) {
	var maxInt = 0;
	for(var peak of peaks) {
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
        var options = container.data("options");
	var me = parseFloat($(getElementSelector(container, elementIds.massError)).val());
	var unit = getMassErrorUnit(container);
	if(me !== options.massError) {
	    options.massError = me;
	    container.data("massErrorChanged", true);
	}
	else if(options.massErrorUnit !== unit)
	{
	    options.massErrorUnit = unit;
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

        var precursorSelectionWin = [];
        if(options.selWinLow && options.selWinHigh)
        {
            precursorSelectionWin = [{ color:"#ccd8e2", xaxis:{from:options.selWinLow, to:options.selWinHigh}, }];
        }
	var ms1plotOptions = {
	    series: { peaks: {show: true, shadowSize: 0}, shadowSize: 0},
	    grid: { show: true,
		    hoverable: true,
		    autoHighlight: true,
		    clickable: true,
                    markings: precursorSelectionWin,
		    borderWidth: 1,
		    labelMargin: 1},
	    selection: { mode: "xy", color: "#F0E68C" },
	    xaxis: { tickLength: 2, tickColor: "#000" },
	    yaxis: { tickLength: 0, tickColor: "#fff", ticks: [], labelWidth: labelWidth }
	};

	if(ms1zoomRange) {
	    ms1plotOptions.xaxis.min = ms1zoomRange.xaxis.from;
	    ms1plotOptions.xaxis.max = ms1zoomRange.xaxis.to;
	    ms1plotOptions.yaxis.min = 0; // ms1zoomRange.yaxis.from;
	    ms1plotOptions.yaxis.max = ms1zoomRange.yaxis.to;
	}

	var placeholder = $(getElementSelector(container, elementIds.msPlot));
	var ms1plot = $.plot(placeholder, data, ms1plotOptions);
        container.data("ms1plot", ms1plot);


	// Mark the precursor peak with a green triangle.
	if(options.precursorMz) {

            var o = ms1plot.pointOffset({ x: options.precursorMz, y: options.precursorIntensity});
            var ctx = ms1plot.getCanvas().getContext("2d");
            ctx.beginPath();
            ctx.moveTo(o.left-10, o.top-5);
            ctx.lineTo(o.left-10, o.top + 5);
            ctx.lineTo(o.left-10 + 10, o.top);
            ctx.lineTo(o.left-10, o.top-5);
            ctx.fillStyle = "#008800";
            ctx.fill();
            placeholder.append('<div style="position:absolute;left:' + (o.left + 4) + 'px;top:' + (o.top-4) + 'px;color:#000;font-size:smaller">'+options.precursorMz.toFixed(2)+'</div>');

	}

	// mark the scan number if we have it
	o = ms1plot.getPlotOffset();
	if(options.ms1scanLabel) {
	    placeholder.append('<div style="position:absolute;left:' + (o.left + 4) + 'px;top:' + (o.top+4) + 'px;color:#666;font-size:smaller">MS1 scan: '+options.ms1scanLabel+'</div>');
	}

	// zoom out icon on plot right hand corner if we are not already zoomed in to the precursor.
	if(container.data("ms1zoomRange")) {
	    placeholder.append('<div id="'+getElementId(container, elementIds.ms1plot_zoom_out)+'" class="zoom_out_link"  style="position:absolute; left:'
			       + (o.left + ms1plot.width() - 40) + 'px;top:' + (o.top+4) + 'px;"></div>');

	    $(getElementSelector(container, elementIds.ms1plot_zoom_out)).click( function() {
		container.data("ms1zoomRange", null);
		createMs1Plot(container);
	    });
	}
	else {
	    placeholder.append('<div id="'+getElementId(container, elementIds.ms1plot_zoom_in)+'" class="zoom_in_link"  style="position:absolute; left:'
			       + (o.left + ms1plot.width() - 20) + 'px;top:' + (o.top+4) + 'px;"></div>');
	    $(getElementSelector(container, elementIds.ms1plot_zoom_in)).click( function() {
		var ranges = {};
		ranges.yaxis = {};
		ranges.xaxis = {};
		ranges.yaxis.from = 0.0;
		ranges.yaxis.to = options.maxIntensityInMs1ZoomRange;

                ranges.xaxis.from = options.precursorMz - 5.0;
                ranges.xaxis.to = options.precursorMz + 5.0;

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
	if(!container.data("zoomRange"))
        {
            plot = $.plot($(getElementSelector(container, elementIds.msmsplot)), datasets,  container.data("plotOptions"));
        }
	else {
            var zoomRange = container.data("zoomRange");
            var selectOpts = {};
	    if($(getElementSelector(container, elementIds.zoom_x)).is(":checked"))
		selectOpts.xaxis = { min: zoomRange.xaxis.from, max: zoomRange.xaxis.to };
	    if($(getElementSelector(container, elementIds.zoom_y)).is(":checked"))
		selectOpts.yaxis = { min: 0, max: zoomRange.yaxis.to };

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

        // Draw the peak mass error plot
        plotPeakMassErrorPlot(container, datasets);
        if(container.data("options").showMassErrorPlot === false)
        {
            $(getElementSelector(container, elementIds.massErrorPlot)).hide();
        }
    }

    function getPlotXRange(options) {

        var xmin = options.minDisplayMz;
        var xmax = options.maxDisplayMz;
        var xpadding = (xmax - xmin) * 0.025;
        // console.log("x-axis padding: "+xpadding);
        return {xmin:xmin - xpadding, xmax:xmax + xpadding};
    }

    function plotPeakMassErrorPlot(container, datasets) {
        var data = [];
        var options = container.data("options");

        var ppmError = options.massErrorPlotDefaultUnit === massErrorTypePpm;

        var minMassError = 0;
        var maxMassError = 0;
        for (var i = 0; i < datasets.length; i += 1) {
            var series = datasets[i];
            var seriesType = series.labelType;
            if (seriesType && seriesType === 'ion') {
                var seriesData = series.data;
                var s_data = [];
                for (var j = 0; j < seriesData.length; j += 1) {
                    var observedMz = seriesData[j][0];
                    var theoreticalMz = seriesData[j][2];
                    var massError = theoreticalMz - observedMz;
                    if (ppmError) {
                        massError = ((massError) / observedMz) * 1000000;
                    }
                    minMassError = Math.min(minMassError, massError);
                    maxMassError = Math.max(maxMassError, massError);
                    s_data.push([seriesData[j][0], massError]);
                }
                data.push({data:s_data, color:series.color, labelType:'none'});
            }
        }

        var placeholder = $(getElementSelector(container, elementIds.massErrorPlot));

        var __xrange = getPlotXRange(options);
        var zoomRange = container.data("zoomRange");
        if(zoomRange)
        {
            // Sync zooming with the MS/MS plot.
            __xrange.xmin = zoomRange.xaxis.from;
            __xrange.xmax = zoomRange.xaxis.to;
        }

        var ypadding = Math.abs(maxMassError - minMassError) * 0.025;

        // the MS/MS plot should have been created by now.  This is a hack to get the plots aligned.
        // We will set the y-axis labelWidth to this value.
        var labelWidth = container.data("plot").getAxes().yaxis.labelWidth;

        if(data.length === 0)
        {
            // Add dummy data to show an empty plot.
            data.push({data:[__xrange.xmin, 0], color:'#000000', labelType:'none'});
        }

        var massErrorPlotOptions = {
            series:{data:data, points:{show:true, fill:true, radius:1}, shadowSize:0},
            grid:{ show:true,
                   hoverable:true,
                   autoHighlight:false,
                   clickable:false,
                   borderWidth:1,
                   labelMargin:1,
                   markings:[
                       {yaxis:{from:0, to:0}, color:"#555555", lineWidth:0.5}
                   ]  // draw a horizontal line at y=0
		 },
            selection:{ mode:"xy", color:"#F0E68C" },
            xaxis:{ tickLength:3, tickColor:"#000",
                    min:__xrange.xmin,
                    max:__xrange.xmax},
            yaxis:{ tickLength:0, tickColor:"#fff",
                    min:minMassError - ypadding,
                    max:maxMassError + ypadding,
                    labelWidth:labelWidth }
        };


        // TOOLTIPS
        $(getElementSelector(container, elementIds.massErrorPlot)).bind("plothover", function (event, pos, item) {
            displayTooltip(item, container, options, "m/z", "error");
        });

        var massErrorPlot = $.plot(placeholder, data, massErrorPlotOptions);

        // Display clickable mass error unit.
        var o = massErrorPlot.getPlotOffset();
        placeholder.append('<div id="' + getElementId(container, elementIds.massErrorPlot_unit) + '" class="link"  '
			   + 'style="position:absolute; left:'
			   + (o.left + 5) + 'px;top:' + (o.top + 4) + 'px;'
			   + 'background-color:yellow; font-style:italic">'
			   + options.massErrorPlotDefaultUnit + '</div>');
        // Toggle mass error unit on click.
        $(getElementSelector(container, elementIds.massErrorPlot_unit)).click(function () {
            var unit = $(this).text();

            if (unit === massErrorTypeTh) {
                $(this).text(massErrorTypePpm);
                options.massErrorPlotDefaultUnit = massErrorTypePpm;
            }
            else if (unit === massErrorTypePpm) {
                $(this).text(massErrorTypeTh);
                options.massErrorPlotDefaultUnit = massErrorTypeTh;
            }
            plotPeakMassErrorPlot(container, datasets);
        });
    }

    function displayTooltip(item, container, options, tooltip_xlabel, tooltip_ylabel) {

        if ($(getElementSelector(container, elementIds.enableTooltip) + ":checked").length > 0) {
            if (item) {
                if (container.data("previousPoint") != item.datapoint) {
                    container.data("previousPoint", item.datapoint);

                    $(getElementSelector(container, elementIds.msmstooltip)).remove();
                    var x = item.datapoint[0].toFixed(2),
                    y = item.datapoint[1].toFixed(2);

                    showTooltip(container, item.pageX, item.pageY,
				tooltip_xlabel + ": " + x + "<br>" + tooltip_ylabel + ": " + y, options);
                }
            }
            else {
                $(getElementSelector(container, elementIds.msmstooltip)).remove();
                container.data("previousPoint", null);
            }
        }
    }

    // -----------------------------------------------
    // SET UP INTERACTIVE ACTIONS FOR MS/MS PLOT
    // -----------------------------------------------
    function setupInteractions (container, options) {

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
	    plotAccordingToChoices(container);
	});

	// TOOLTIPS
	$(getElementSelector(container, elementIds.msmsplot)).bind("plothover", function (event, pos, item) {
	    displayTooltip(item, container, options, "m/z", "intensity");
	});
	$(getElementSelector(container, elementIds.enableTooltip)).click(function() {
	    $(getElementSelector(container, elementIds.msmstooltip)).remove();
	});

        // PLOT MASS ERROR CHECKBOX
        $(getElementSelector(container, elementIds.massErrorPlot_option)).click(function() {
	    var plotDiv = $(getElementSelector(container, elementIds.massErrorPlot));
	    if($(this).is(':checked'))
	    {
                plotDiv.show();
	    }
	    else
	    {
                plotDiv.hide();
	    }
        });


	// SHOW / HIDE ION SERIES; UPDATE ON MASS TYPE CHANGE;
	// PEAK ASSIGNMENT TYPE CHANGED; PEAK LABEL TYPE CHANGED
	var ionChoiceContainer = $(getElementSelector(container, elementIds.ion_choice));
	ionChoiceContainer.find("input").click(function() {
	    plotAccordingToChoices(container);
        });

        $(getElementSelector(container, elementIds.immoniumIons)).click(function() {
	    plotAccordingToChoices(container);
        });

        $(getElementSelector(container, elementIds.reporterIons)).click(function() {
	    plotAccordingToChoices(container);
        });

        $(getElementSelector(container, elementIds.labelPrecursor)).click(function() {
	    plotAccordingToChoices(container);
	});

	// Plot neutral loss options
	var neutralLossContainer = $(getElementSelector(container, elementIds.nl_choice));
	neutralLossContainer.find("input").click(function() {
	    container.data("selectedNeutralLossChanged", true);
	    var selectedNeutralLosses = getNeutralLosses(container);
	    container.data("options").peptide.recalculateLossOptions(selectedNeutralLosses, container.data("options").maxNeutralLossCount);
	    plotAccordingToChoices(container);
	});

        // Mass type options
	container.find("input[name='"+getRadioName(container, "massTypeOpt")+"']").click(function() {
	    container.data("massTypeChanged", true);
	    plotAccordingToChoices(container);
	});

        // Peak detect checkbox
        $(getElementSelector(container, elementIds.peakDetect)).click(function() {
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
	makeIonTableMovable(container, options);

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

    function showTooltip(container, x, y, contents, options) {
	var tooltipCSS = {
            position: 'absolute',
            display: 'none',
            top: y + 5,
            left: x + 5,
            border: '1px solid #fdd',
            padding: '2px',
            'background-color': '#F0E68C',
            opacity: 0.80 };

	if ( options.tooltipZIndex !== undefined && options.tooltipZIndex !== null ) {
	    tooltipCSS["z-index"] = options.tooltipZIndex;
	}

        $('<div id="'+getElementId(container, elementIds.msmstooltip)+'">' + contents + '</div>')
	    .css( tooltipCSS ).appendTo("body").fadeIn(200);
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
                $(getElementSelector(container, elementIds.massErrorPlot)).css({width: width});

		plotAccordingToChoices(container);
		if(options.ms1peaks && options.ms1peaks.length > 0) {
		    $(getElementSelector(container, elementIds.msPlot)).css({width: width});
		    createMs1Plot(container);
		}
		$(getElementSelector(container, elementIds.slider_width_val)).text(width);
		if ( options.sizeChangeCallbackFunction ) {
		    options.sizeChangeCallbackFunction();
		}
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
		if ( options.sizeChangeCallbackFunction ) {
		    options.sizeChangeCallbackFunction();
		}
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
        //added two peak lists into data
        var peaks1 = {data: options.peaks,  color: "#bbbbbb", labelType: 'none', yaxis:1};
        var peaks2 = {data: options.peaks2, color: "#bbbbbb", labelType: 'none', yaxis:2, label: options.ms2peaks2Label};
        var data = [];
        data.push(peaks1);
        data.push(peaks2);


	// add the annotated peaks
	var seriesMatches = getSeriesMatches(container, selectedIonTypes);
	for(var i = 0; i < seriesMatches.length; i += 1) {
	    data.push(seriesMatches[i]);
	}

        // add immonium ions
        if(labelImmoniumIons(container))
        {
            // Recalculate if the following have changed: mass error, peak assignment type
            if(container.data("massErrorChanged") || container.data("peakAssignmentTypeChanged"))
            {
                calculateImmoniumIons(container);
            }
            data.push(container.data("immoniumIons"));
        }

        // add precursor peak
        if(labelPrecursorPeak(container))
        {
            // Recalculate if the following have changed: mass error, peak assignment type, selected neutral loss
            // Changing mass type only chanages the fragmentMassType not the precursorMassType, so the labeled precursor peaks will not change.
            if(container.data("massErrorChanged") || container.data("peakAssignmentTypeChanged") || container.data("selectedNeutralLossChanged"))
            {
                calculatePrecursorPeak(container);
            }
            if (container.data("precursorPeak"))
            {
                data.push(container.data("precursorPeak"));
            }
        }

        // add reporter ions
        if(labelReporterIons(container))
        {
            // Recalculate if the following have changed: mass error, peak assignment type
            if(container.data("massErrorChanged") || container.data("peakAssignmentTypeChanged"))
            {
                calculateReporterIons(container);
            }
            var reporterSeries = container.data("reporterSeries");
            for(var i = 0; i < reporterSeries.length; i += 1)
            {
                var matches = reporterSeries[i].matches;
                if(matches)  data.push(matches);
            }
        }

        // add user reporter ions
        if(container.data("userReporterSeries")) {
            data.push(container.data("userReporterSeries"));
        }

	// add any user specified extra peaks
	for(var i = 0; i < options.extraPeakSeries.length; i += 1) {
	    data.push(options.extraPeakSeries[i]);
	}
	return data;
    }

    function labelPrecursorPeak(container)
    {
        return $(getElementSelector(container, elementIds.labelPrecursor)).is(":checked")
    }

    function labelImmoniumIons(container)
    {
        return $(getElementSelector(container, elementIds.immoniumIons)).is(":checked")
    }

    function labelReporterIons(container)
    {
        return $(getElementSelector(container, elementIds.reporterIons)).is(":checked");
    }

    function calculateImmoniumIons(container)
    {
        var options = container.data("options");
        var antic = 0;

        var peaks = options.peaks;

        // immonium ions (70 P, 72 V, 86 I/L, 110 H, 120 F, 136 Y, 159 W)
        var immoniumIonTypes = [{mass:70.0657, aa:'P'},
                                {mass:102.0555, aa:'E'},
                                {mass:72.0813, aa:'V'},
                                {mass:86.0970, aa:'I/L'},
                                {mass:110.0718, aa:'H'},
                                {mass:120.0813, aa:'F'},
                                {mass:101.0715, aa:'Q'},
                                {mass:136.0762, aa:'Y'},
                                {mass:159.0922, aa:'W'}];

        var immoniumIonMatches = [];
        var labels = [];
        for(var i = 0; i < immoniumIonTypes.length; i += 1)
        {
            var ion = immoniumIonTypes[i];
            var match = getMatchingPeakForMz(container, peaks, ion.mass);
            if(match.bestPeak)
            {
                immoniumIonMatches.push([match.bestPeak[0], match.bestPeak[1]]);
                labels.push(ion.aa + '-' + match.bestPeak[0].toFixed(1));
		antic += match.bestPeak[1];
            }
        }
        container.data("anticImmonium", antic);
        container.data("immoniumIons", {data: immoniumIonMatches, labels: labels, color: "#008000"});
    }

    function calculatePrecursorPeak(container) {
        var options = container.data("options");
        var antic = 0;

        if(options.labelPrecursorPeak && options.theoreticalMz) {
            var precursorMzMatches = [];
            var labels = [];

            var peaks = options.peaks;
            var precursorMz = options.theoreticalMz;
            var charge = options.charge ? options.charge : 1;
            var label = 'M';

            // Label all possible precursor charge states
            for (var i = 1; i <= charge; i += 1) {
                label += "+";
                var pmz = (((precursorMz - MASS_PROTON) * charge) / i) + MASS_PROTON;
                //console.log("for charge:"+i+" -- m/z:"+pmz);

                var match = getMatchingPeakForMz(container, peaks, pmz);
                if (match.bestPeak) {
                    precursorMzMatches.push([match.bestPeak[0], match.bestPeak[1]]);
                    labels.push(label);
                    antic += match.bestPeak[1];
                }

                var neutralLosses = getNeutralLosses(container);
                for (var lossKey in neutralLosses) {
                    var loss = neutralLosses[lossKey];
                    // console.log("Neutral loss:"+lossKey+" -- Label:"+loss.label());

                    match = getMatchingPeakForMz(container, peaks, pmz - (loss.monoLossMass / charge));
                    if (match.bestPeak) {
                        precursorMzMatches.push([match.bestPeak[0], match.bestPeak[1]]);
                        labels.push(label + " " + loss.label());
                        options.antic += match.bestPeak[1];
                    }
                }
                if (precursorMzMatches.length > 0)
                    container.data("precursorPeak", {data: precursorMzMatches, labels: labels, color: "#ffd700"});
            }
        }
        container.data("anticPrecursor", antic);
    }

    function calculateReporterIons(container)
    {
        container.data("anticReporter", 0);

        var itraqWholeLabel = 145.1069;
        var itraqIons = [113.107325, 114.11068, 115.107715, 116.111069, 117.114424, 118.111459, 119.114814, 121.121524];

        var tmtWholeLabel = 230.1702;
        var tmtIons = [126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.137790, 130.134825, 130.141145, 131.138180, 131.144500, 132.141535, 132.147855, 133.14489, 133.15121, 134.148245]; // tmt-16

        var reporterSeries = [];
        reporterSeries.push({color: "#2f4f4f", ions: itraqIons, wholeLabel: itraqWholeLabel});  // DarkSlateBlue
        reporterSeries.push({color: "#556b2f", ions: tmtIons, wholeLabel: tmtWholeLabel});  // DarkOliveGreen

        for(var i = 0; i < reporterSeries.length; i += 1)
        {
            var series = reporterSeries[i];
            var matches = calculateReporters(series, container);
            reporterSeries[i].matches = matches;
        }

        container.data("reporterSeries", reporterSeries);
    }

    /**
     * Calculate user-supplied reporter ion matches.
     *
     * @param container
     */
    function calculateUserReporters(options, container) {

        const ionMzArray = options.userReporterIons;
        let antic = 0;

        const peaks = options.peaks;
        let matches = [];
        let labels = [];

        for(let i = 0; i < ionMzArray.length; i++)
        {
            const match = getMatchingPeakForMz(container, peaks, ionMzArray[i]);
            if(match.bestPeak)
            {
                matches.push([match.bestPeak[0], match.bestPeak[1]]);
                labels.push(match.bestPeak[0].toFixed(2));
                antic += match.bestPeak[1];
            }
        }

        if(matches.length > 0)
        {
            container.data("userReporterSeries",  {data: matches, labels: labels, color:"#2f4f4f"} ); // DarkSlateBlue
            container.data("anticUserReporterIons", antic);
        }
    }


    function calculateReporters(seriesInfo, container)
    {
        var ionMzArray = seriesInfo.ions;

        var options = container.data("options");
        var peaks = options.peaks;
        var matches = [];
        var antic = container.data("anticReporter");
        for(var i = 0; i < ionMzArray.length; i += 1)
        {
            var match = getMatchingPeakForMz(container, peaks, ionMzArray[i]);
            if(match.bestPeak)
            {
                matches.push([match.bestPeak[0], match.bestPeak[1]]);
		antic += match.bestPeak[1];
            }
        }
        var labels = [];
        var maxIntensity = 0;
        for(var i = 0; i < matches.length; i += 1)
        {
            maxIntensity = Math.max(maxIntensity, matches[i][1]);
        }

        for(var i = 0; i < matches.length; i += 1)
        {
            var intensity = matches[i][1];
            var rank = Math.round(((intensity/maxIntensity) * 100.0));
            var mzRounded = matches[i][0].toFixed(1);
            labels.push(mzRounded + " (" + rank + "%)");
        }
        // Get a match for the whole label.
        var match = getMatchingPeakForMz(container, peaks, seriesInfo.wholeLabel);
        if(match.bestPeak)
        {
            matches.push([match.bestPeak[0], match.bestPeak[1]]);
            labels.push('nterm');
	    antic += match.bestPeak[1];
        }

        container.data("anticReporter", antic);
        if(matches.length > 0)
        {
            return {data: matches, labels: labels, color:seriesInfo.color};
        }
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
            {
		ntermIons.push(sion);
	    }
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
    // CALCULATE THEORETICAL MASSES FOR THE SELECTED ION SERIES
    // ---------------------------------------------------------
    function calculateTheoreticalSeries(container, selectedIons) {

        if(container.data("massTypeChanged"))
        {
            // Clear out the theoretical ion series if the selected mass type changed.
            container.data("ionSeries", {a: [], b: [], c: [], x: [], y: [], z: []});
        }
	if(selectedIons) {

	    var todoIonSeries = [];
	    var todoIonSeriesData = [];
            var ionSeries = container.data("ionSeries");
	    for(var i = 0; i < selectedIons.length; i += 1) {
		var sion = selectedIons[i];
		if(sion.type == "a") {
		    if(ionSeries.a[sion.charge])	continue; // already calculated
		    else {
			todoIonSeries.push(sion);
			ionSeries.a[sion.charge] = [];
			todoIonSeriesData.push(ionSeries.a[sion.charge]);
		    }
		}
		if(sion.type == "b") {
		    if(ionSeries.b[sion.charge])	continue; // already calculated
		    else {
			todoIonSeries.push(sion);
			ionSeries.b[sion.charge] = [];
			todoIonSeriesData.push(ionSeries.b[sion.charge]);
		    }
		}
		if(sion.type == "c") {
		    if(ionSeries.c[sion.charge])	continue; // already calculated
		    else {
			todoIonSeries.push(sion);
			ionSeries.c[sion.charge] = [];
			todoIonSeriesData.push(ionSeries.c[sion.charge]);
		    }
		}
		if(sion.type == "x") {
		    if(ionSeries.x[sion.charge])	continue; // already calculated
		    else {
			todoIonSeries.push(sion);
			ionSeries.x[sion.charge] = [];
			todoIonSeriesData.push(ionSeries.x[sion.charge]);
		    }
		}
		if(sion.type == "y") {
		    if(ionSeries.y[sion.charge])	continue; // already calculated
		    else {
			todoIonSeries.push(sion);
			ionSeries.y[sion.charge] = [];
			todoIonSeriesData.push(ionSeries.y[sion.charge]);
		    }
		}
		if(sion.type == "z") {
		    if(ionSeries.z[sion.charge])	continue; // already calculated
		    else {
			todoIonSeries.push(sion);
			ionSeries.z[sion.charge] = [];
			todoIonSeriesData.push(ionSeries.z[sion.charge]);
		    }
		}
	    }

	    if(container.data("options").sequence) {

                var sequence = container.data("options").sequence
		var massType = getMassType(container);

		for(var i = 1; i < sequence.length; i += 1) {

		    for(var j = 0; j < todoIonSeries.length; j += 1) {
			var tion = todoIonSeries[j];
			var ionSeriesData = todoIonSeriesData[j];

                        var ion = Ion.getSeriesIon(tion, container.data("options").peptide, i, massType);
                        // Put the ion masses in increasing value of m/z, For c-term ions the array will have to be
                        // populated backwards.
			if(tion.term == "n")
                            // Add to end of array
			    ionSeriesData.push(ion);
			else if(tion.term == "c")
                            // Add to beginning of array
			    ionSeriesData.unshift(ion);
		    }
		}
	    }
	}
    }

    function getMassType(container)
    {
        return container.find("input[name='"+getRadioName(container, "massTypeOpt")+"']:checked").val();
    }

    function getMassErrorUnit(container)
    {
        return container.find("input[name='"+getRadioName(container, "massErrorUnit")+"']:checked").val();
    }

    function getPeakAssignmentType(container)
    {
        return container.find("input[name='"+getRadioName(container, "peakAssignOpt")+"']:checked").val();
    }

    // -----------------------------------------------
    // MATCH THEORETICAL MASSES WITH PEAKS IN THE SCAN
    // -----------------------------------------------
    function recalculateFragmentMatches(container) {
	return (container.data("massErrorChanged") ||
		container.data("massTypeChanged") ||
		container.data("peakAssignmentTypeChanged") ||
		container.data("selectedNeutralLossChanged"));
    }

    function clearIonSeries(container)
    {
        container.data("ionSeriesMatch", {a: [], b: [], c: [], x: [], y: [], z: []});
        container.data("ionSeriesLabels", {a: [], b: [], c: [], x: [], y: [], z: []});
        container.data("ionSeriesAntic", {a: [], b: [], c: [], x: [], y: [], z: []});
    }

    function getSeriesMatches(container, selectedIonTypes) {

	var dataSeries = [];

        if(recalculateFragmentMatches(container))
        {
            // If mass type, mass error, peak assignment or selected neutral losses have changed
            // all peak series should be recalculated.
            clearIonSeries(container);
        }

	var peakAssignmentType = getPeakAssignmentType(container);
	var peakLabelType = container.find("input[name='"+getRadioName(container, "peakLabelOpt")+"']:checked").val();
        var massType = getMassType(container);
        var massErrorUnit = getMassErrorUnit(container);

        var ionSeriesMatch = container.data("ionSeriesMatch");
        var ionSeries = container.data("ionSeries");
        var ionSeriesLabels = container.data("ionSeriesLabels");
        var ionSeriesAntic = container.data("ionSeriesAntic");
        var massError = container.data("options").massError;
        var peaks = getPeaks(container);



	for(var j = 0; j < selectedIonTypes.length; j += 1) {

	    var ion = selectedIonTypes[j];

	    if(ion.type == "a") {
		if(!ionSeriesMatch.a[ion.charge]) { // re-calculate only if matching peaks for this series have not been calculated already
		    // calculated matching peaks
		    var adata = calculateMatchingPeaks(container, ionSeries.a[ion.charge], peaks, massError, massErrorUnit, peakAssignmentType, massType);
		    if(adata && adata.length > 0) {
			ionSeriesMatch.a[ion.charge] = adata[0];
			ionSeriesLabels.a[ion.charge] = adata[1];
			ionSeriesAntic.a[ion.charge] = adata[2];
		    }
		}
		dataSeries.push({data: ionSeriesMatch.a[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.a[ion.charge]});
	    }

	    if(ion.type == "b") {
		if(!ionSeriesMatch.b[ion.charge]) { // re-calculate only if matching peaks for this series have not been calculated already
		    // calculated matching peaks
		    var bdata = calculateMatchingPeaks(container, ionSeries.b[ion.charge], peaks, massError, massErrorUnit, peakAssignmentType, massType);
		    if(bdata && bdata.length > 0) {
			ionSeriesMatch.b[ion.charge] = bdata[0];
			ionSeriesLabels.b[ion.charge] = bdata[1];
                        ionSeriesAntic.b[ion.charge] = bdata[2];
		    }
		}
		dataSeries.push({data: ionSeriesMatch.b[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.b[ion.charge]});
	    }

	    if(ion.type == "c") {
		if(!ionSeriesMatch.c[ion.charge]) { // re-calculate only if matching peaks for this series have not been calculated already
		    // calculated matching peaks
		    var cdata = calculateMatchingPeaks(container, ionSeries.c[ion.charge], peaks, massError, massErrorUnit, peakAssignmentType, massType);
		    if(cdata && cdata.length > 0) {
			ionSeriesMatch.c[ion.charge] = cdata[0];
			ionSeriesLabels.c[ion.charge] = cdata[1];
                        ionSeriesAntic.c[ion.charge] = cdata[2];
		    }
		}
		dataSeries.push({data: ionSeriesMatch.c[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.c[ion.charge]});
	    }

	    if(ion.type == "x") {
		if(!ionSeriesMatch.x[ion.charge]) { // re-calculate only if matching peaks for this series have not been calculated already
		    // calculated matching peaks
		    var xdata = calculateMatchingPeaks(container, ionSeries.x[ion.charge], peaks, massError, massErrorUnit, peakAssignmentType, massType);
		    if(xdata && xdata.length > 0) {
			ionSeriesMatch.x[ion.charge] = xdata[0];
			ionSeriesLabels.x[ion.charge] = xdata[1];
                        ionSeriesAntic.x[ion.charge] = xdata[2];
		    }
		}
		dataSeries.push({data: ionSeriesMatch.x[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.x[ion.charge]});
	    }

	    if(ion.type == "y") {
		if(!ionSeriesMatch.y[ion.charge]) { // re-calculate only if matching peaks for this series have not been calculated already
		    // calculated matching peaks
		    var ydata = calculateMatchingPeaks(container, ionSeries.y[ion.charge], peaks, massError, massErrorUnit, peakAssignmentType, massType);
		    if(ydata && ydata.length > 0) {
			ionSeriesMatch.y[ion.charge] = ydata[0];
			ionSeriesLabels.y[ion.charge] = ydata[1];
                        ionSeriesAntic.y[ion.charge] = ydata[2];
		    }
		}
		dataSeries.push({data: ionSeriesMatch.y[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.y[ion.charge]});
	    }

	    if(ion.type == "z") {
		if(!ionSeriesMatch.z[ion.charge]) { // re-calculate only if matching peaks for this series have not been calculated already
		    // calculated matching peaks
		    var zdata = calculateMatchingPeaks(container, ionSeries.z[ion.charge], peaks, massError, massErrorUnit, peakAssignmentType, massType);
		    if(zdata && zdata.length > 0) {
			ionSeriesMatch.z[ion.charge] = zdata[0];
			ionSeriesLabels.z[ion.charge] = zdata[1];
                        ionSeriesAntic.z[ion.charge] = zdata[2];
		    }
		}
		dataSeries.push({data: ionSeriesMatch.z[ion.charge], color: ion.color, labelType: peakLabelType, labels: ionSeriesLabels.z[ion.charge]});
	    }
	}
	return dataSeries;
    }

    function getNeutralLosses(container) {
        var neutralLosses = [];
        $(getElementSelector(container, elementIds.nl_choice)).find("input:checked").each(function() {
            var lossLabel = $(this).val();
            var loss = container.data("options").peptide.getLossForLabel(lossLabel);
            neutralLosses.push(loss);
        });
        return neutralLosses;
    }

    function getLabel(sion, neutralLosses) {
        var label = sion.label;
        if(neutralLosses) {
            label += neutralLosses.getLabel();
        }
        return label;
    }

    function ionMz(sion, neutralLosses, massType) {
        var ionmz;
        if(!neutralLosses)
            ionmz = sion.mz;
        else {
            ionmz = Ion.getIonMzWithLoss(sion, neutralLosses, massType);
        }
        return ionmz;
    }

    function calculateMatchingPeaks(container, ionSeries, allPeaks, massTolerance, massErrorUnit, peakAssignmentType, massType) {

        // console.log("calculating matching peaks");
	var peakIndex = 0;

	var matchData = [];
	matchData[0] = []; // peaks
	matchData[1] = []; // labels -- ions;
        matchData[2] = 0;  // Annotated ion current for this series.
	var peptide = container.data("options").peptide;

	for(var i = 0; i < ionSeries.length; i += 1) {

	    var sion = ionSeries[i];

	    // get match for water and / or ammonia loss.
            var minIndex = Number.MAX_VALUE;
            var neutralLossOptions = peptide.getPotentialLosses(sion);

	    var index = getMatchForIon(sion, matchData, allPeaks, peakIndex, massTolerance, massErrorUnit, peakAssignmentType, null, massType);
            minIndex = Math.min(minIndex, index);

            for(var n = 1; n < neutralLossOptions.length; n += 1)
            {
                var loss_options_with_n_losses = neutralLossOptions[n];
                for(var k = 0; k < loss_options_with_n_losses.lossCombinationCount(); k += 1)
                {
                    var lossCombination = loss_options_with_n_losses.getLossCombination(k);
		    var index = getMatchForIon(sion, matchData, allPeaks, peakIndex, massTolerance, massErrorUnit, peakAssignmentType, lossCombination, massType);
                    minIndex = Math.min(minIndex, index);
                }
            }

            peakIndex = minIndex;
	}

	return matchData;
    }

    // sion -- theoretical ion
    // matchData -- array to which we will add a peak if there is a match
    // allPeaks -- array with all the scan peaks
    // peakIndex -- current index in peaks array
    // Returns the index of the matching peak, if one is found
    function getMatchForIon(sion, matchData, allPeaks, peakIndex, massTolerance, massErrorUnit, peakAssignmentType, neutralLosses, massType) {

	if(!neutralLosses)
	    sion.match = false; // reset;
	var ionmz = ionMz(sion, neutralLosses, massType);
        var peakLabel = getLabel(sion, neutralLosses);

	var __ret = getMatchingPeak(peakIndex, allPeaks, ionmz, massTolerance, massErrorUnit, peakAssignmentType);

        peakIndex = __ret.peakIndex;
        var bestPeak = __ret.bestPeak;

        // if we found a matching peak for the current ion, save it
        if(bestPeak) {
            // console.log("found match "+sion.label+", "+ionmz+";  peak: "+bestPeak[0] + "; theoreticalMz: " + __ret.theoreticalMz);
            matchData[0].push([bestPeak[0], bestPeak[1], __ret.theoreticalMz]);
            matchData[1].push(peakLabel);
            matchData[2] = matchData[2] + bestPeak[1];
            if(!neutralLosses) {
                sion.match = true;
            }
        }

	return peakIndex;
    }

    function getMatchingPeakForMz(container, allPeaks, ionMz)
    {
        var massError = container.data("options").massError;
        var massErrorUnit = getMassErrorUnit(container);
        var peakAssignmentType = getPeakAssignmentType(container);
        return getMatchingPeak(0, allPeaks, ionMz, massError, massErrorUnit, peakAssignmentType);
    }

    function getMatchingPeak(peakIndex, allPeaks, ionmz, massTolerance, toleranceUnit, peakAssignmentType) {

        var bestDistance;
        var bestPeak;
        var tolerantPeakMin = ionmz - massTolerance;
        var tolerantPeakMax = ionmz + massTolerance;

        for (var j = peakIndex; j < allPeaks.length; j += 1) {

            var peak = allPeaks[j];

            if(toleranceUnit === massErrorTypePpm)
            {
                var tolerance = (massTolerance * peak[0])/1000000;
                tolerantPeakMin = ionmz - tolerance;
                tolerantPeakMax = ionmz + tolerance;
            }
            // peak is before the current ion we are looking at
            if (peak[0] < tolerantPeakMin)
                continue;

            // peak is beyond the current ion we are looking at
            if (peak[0] > tolerantPeakMax) {
                peakIndex = j;
                break;
            }

            // peak is within +/- massTolerance of the current ion we are looking at

            // if this is the first peak in the range
            if (!bestPeak) {
                //console.log("found a peak in range, "+peak.mz);
                bestPeak = peak;
                bestDistance = Math.abs(ionmz - peak[0]);
                continue;
            }

            // if peak assignment method is Most Intense
            if (peakAssignmentType == "intense") {
                if (peak[1] > bestPeak[1]) {
                    bestPeak = peak;
                    continue;
                }
            }

            // if peak assignment method is Closest Peak
            if (peakAssignmentType == "close") {
                var dist = Math.abs(ionmz - peak[0]);
                if (!bestDistance || dist < bestDistance) {
                    bestPeak = peak;
                    bestDistance = dist;
                }
            }
        }

        return {peakIndex:peakIndex, bestPeak:bestPeak, theoreticalMz:ionmz};
    }


    function getPeaks(container)
    {
        var options = container.data("options");

        if($(getElementSelector(container, elementIds.peakDetect)).is(":checked"))
        {
            if(options.sparsePeaks == null) {
                doPeakDetection(container);
            }
            return options.sparsePeaks;
        }
        else
        {
            return options.peaks;
        }
    }

    function doPeakDetection(container) {

        // console.log("calculating sparse peaks");

        var peaks = container.data("options").peaks;
        var sparsePeaks = [];

        var intensities = [];
        for(var i = 0; i < peaks.length; i+= 1)
        {
            intensities.push(peaks[i][1]);
        }
        intensities.sort(function(a,b){return b-a});
        var max_50_intensity = intensities[Math.min(intensities.length - 1, 49)];


        var window = 50.0;
        for(var i = 0; i < peaks.length; i += 1) {

	    var peak = peaks[i];

            var intensity = peak[1];

            // If this is one of the 50 most intense peaks, add it to sparse peaks
            if(intensity >= max_50_intensity)
            {
                sparsePeaks.push(peak);
                continue;
            }

            var mz = peak[0];
            var j = i-1;
            var totalIntensity = intensity;
            var peakCount = 1;
            // sum up the intensities in the +/- 50Da window of this peak
            var maxIntensity = intensity;
            var minIndex = i;
            var maxIndex = i;
            while(j >= 0)
            {
                if(peaks[j][0] < mz - window)
                    break;

                if(peaks[j][1] > maxIntensity)
                {
                    maxIntensity = peaks[j][1];
                }
                totalIntensity += peaks[j][1];
                minIndex = j;
                j -= 1;
                peakCount += 1;
            }
            j = i+1;
            while(j < peaks.length)
            {
                if(peaks[j][0] > mz + window)
                    break;

                if(peaks[j][1] > maxIntensity)
                {
                    maxIntensity = peaks[j][1];
                }
                totalIntensity += peaks[j][1];
                maxIndex = j;
                j += 1;
                peakCount += 1;
            }

            var mean = totalIntensity / peakCount;
	    if(peakCount <= 2 && intensity >= maxIntensity * 0.8)
	    {
		sparsePeaks.push(peak);
	    }
	    else if(peakCount <= 10 && intensity == maxIntensity)
            {
                sparsePeaks.push(peak);
            }

            else
            {
                // calculate the standard deviation
                var sdev = 0;
                for(var k = minIndex; k <= maxIndex; k += 1)
                {
                    sdev += Math.pow((peaks[k][1] - mean), 2);
                }
                sdev = Math.sqrt(sdev / peakCount);

                if(intensity >= mean + 2 * sdev)
                {
                    sparsePeaks.push(peak);
                }
                // console.log(intensity+"  "+mean+"  "+sdev);
            }
	}
        // console.log("Sparse Peak count: "+sparsePeaks.length);
        // console.log("All Peaks count: "+peaks.length);
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

	var parentTable = '<table cellpadding="0" cellspacing="5" class="lorikeet-outer-table"> ';
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

        if(options.showSequenceInfo) {
            // placeholder for sequence, m/z, scan number etc
            parentTable += '<td style="background-color: white; padding:5px; border:1px dotted #cccccc;" valign="bottom" align="center"> ';
            parentTable += '<div id="'+getElementId(container, elementIds.seqinfo)+'" style="width:100%;"></div> ';
            // placeholder for file name, scan number and charge
            parentTable += '<div id="'+getElementId(container, elementIds.fileinfo)+'" style="width:100%;"></div> ';
            parentTable += '</td> ';
        }


        if(options.showIonTable) {
            // placeholder for the ion table
            parentTable += '<td rowspan="'+rowspan+'" valign="top" id="'+getElementId(container, elementIds.ionTableLoc1)+'" > ';
            parentTable += '<div id="'+getElementId(container, elementIds.ionTableDiv)+'">';
            parentTable += '<span id="'+getElementId(container, elementIds.moveIonTable)+'" class="font_small link">[Click]</span> <span class="font_small">to move table</span>';
            // placeholder for annotated ion current value
            parentTable += '<div id="'+getElementId(container, elementIds.anticInfo)+'" style="margin-top:5px;"></div> ';
            // placeholder for modifications
            parentTable += '<div id="'+getElementId(container, elementIds.modInfo)+'" style="margin-top:20px;"></div> ';
            // placeholder for user-supplied reporter ion masses
            parentTable += '<div id="'+getElementId(container, elementIds.userReporterIons)+'" style="margin-top:20px;"></div> ';
            parentTable += '</div> ';
            parentTable += '</td> ';
            parentTable += '</tr> ';
        }


	// placeholders for the ms/ms plot
	parentTable += '<tr> ';
	parentTable += '<td style="background-color: white; padding:5px; border:1px dotted #cccccc;" valign="top" align="center"> ';
	parentTable += '<div id="'+getElementId(container, elementIds.msmsplot)+'" align="bottom" style="width:'+options.width+'px;height:'+options.height+'px;"></div> ';

	// placeholder for viewing options (zoom, plot size etc.)
	parentTable += '<div id="'+getElementId(container, elementIds.viewOptionsDiv)+'" align="top" style="margin-top:15px;"></div> ';

        // placeholder for peak mass error plot
        parentTable += '<div id="'+getElementId(container, elementIds.massErrorPlot)+'" style="width:'+options.width+'px;height:100px;"></div> ';

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
	myTable += '<table id="'+getElementId(container, elementIds.ionTable)+'" cellpadding="2" class="font_small '+elementIds.ionTable+'">';
	myTable +=  "<thead>";
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
	myTable +=   "</tr>";
	myTable +=  "</thead>";

	myTable +=  "<tbody>";

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
		    else if(seriesData[i].mz < options.peaks[0][0] |
			    seriesData[i].mz > options.peaks[options.peaks.length - 1][0]) {
			cls="numCell";
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
                    else if(seriesData[idx].mz < options.peaks[0][0] |
                            seriesData[idx].mz > options.peaks[options.peaks.length - 1][0]) {
			cls="numCell";
                    }
		    myTable +=    "<td class='"+cls+"' "+style+" >" +round(seriesData[idx].mz)+  "</td>";
		}
		else {
		    myTable +=    "<td>" +"&nbsp;"+  "</td>";
		}
	    }

	}
	myTable +=   "</tr>";

	myTable += "</tbody>";
	myTable += "</table>";

	// alert(myTable);
	$(getElementSelector(container, elementIds.ionTable)).remove();
	$(getElementSelector(container, elementIds.ionTableDiv)).prepend(myTable); // add as the first child

	if ( options.sizeChangeCallbackFunction ) {
	    options.sizeChangeCallbackFunction();
	}

        showAnticInfo(container); // Total annotated ion current
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

    function makeIonTableMovable(container, options) {

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

	    if ( options.sizeChangeCallbackFunction ) {
		options.sizeChangeCallbackFunction();
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

	    var neutralMass = 0;

	    if(options.precursorMassType == 'mono')
		neutralMass = options.peptide.getNeutralMassMono();
            else
                neutralMass = options.peptide.getNeutralMassAvg();

	    // console.log(options.precursorMassType + " " + neutralMass);
	    var mz;
	    if(options.charge) {
		mz = Ion.getMz(neutralMass, options.charge);
	    }

            // save the theoretical m/z in the options
            options.theoreticalMz = mz;

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
		fileinfo += ', Exp. m/z: '+options.precursorMz;
            }
	    if(options.charge) {
		fileinfo += ', Charge: '+options.charge;
	    }
	    fileinfo += '</div>';
	}

	$(getElementSelector(container, elementIds.fileinfo)).append(fileinfo);
    }

    //---------------------------------------------------------
    // USER-SPECIFIED REPORTER ION INFO
    //---------------------------------------------------------
    function showUserReporterIonsInfo (container) {

        var options = container.data("options");
        var reporterIonInfo = '';

        if(!options.userReporterIons) {
            return;
        }

        reporterIonInfo += '<div>';
        reporterIonInfo += 'Reporter Ions:';
        for( let i = 0; i < options.userReporterIons.length; i++ ) {
            reporterIonInfo += "<div><b>"+ options.userReporterIons[i]+"</b></div>";
        }
        reporterIonInfo += '</div>';

        $(getElementSelector(container, elementIds.userReporterIons)).append(reporterIonInfo);
    }

    //---------------------------------------------------------
    // MODIFICATION INFO
    //---------------------------------------------------------
    function showModInfo (container) {

        var options = container.data("options");

	var modInfo = '';

	modInfo += '<div>';
	if(options.ntermMod || options.ntermMod > 0) {
	    modInfo += 'Add to N-term: <b>'+round(options.ntermMod)+'</b>';
	}
	if(options.ctermMod || options.ctermMod > 0) {
	    modInfo += 'Add to C-term: <b>'+round(options.ctermMod)+'</b>';
	}
	modInfo += '</div>';

	if(options.staticMods && options.staticMods.length > 0) {
	    var sm_modInfo = '<div style="margin-top:5px;">';
	    sm_modInfo += 'Static Modifications: ';
            var count = 0;
	    for(var i = 0; i < options.staticMods.length; i += 1) {
		var mod = options.staticMods[i];
                if(mod.modMass == 0.0)
                    continue;
		//if(i > 0) modInfo += ', ';
		sm_modInfo += "<div><b>"+mod.aa.code+": "+mod.modMass+"</b></div>";
                count += 1;
	    }
	    sm_modInfo += '</div>';
            if(count > 0)
            {
                modInfo += sm_modInfo;
            }
	}

	if(options.variableMods && options.variableMods.length > 0) {

	    var uniqVarMods = {};
	    for(var i = 0; i < options.variableMods.length; i += 1) {
		var mod = options.variableMods[i];
                var varmods = uniqVarMods[mod.aa.code + ' ' + mod.modMass];
		if(!varmods)
                {
		    varmods = [];
                    uniqVarMods[mod.aa.code + ' ' + mod.modMass] = varmods;
                }
		varmods.push(mod);
	    }

            var keys = [];
            for(var key in uniqVarMods)
            {
                if(uniqVarMods.hasOwnProperty(key))
                {
                    keys.push(key);
                }
            }
            keys.sort();

	    modInfo += '<div style="margin-top:5px;">';
	    modInfo += 'Variable Modifications: ';
            modInfo += "<table class='varModsTable'>";
	    for(var k = 0; k < keys.length; k++) {
		var varmods = uniqVarMods[keys[k]];
                modInfo += "<tr><td><span style='font-weight: bold;'>";
                modInfo += varmods[0].aa.code+": "+varmods[0].modMass;
                modInfo += "</span></td>";
                modInfo += "<td>[";
                for(var i = 0; i < varmods.length; i++)
                {
                    if(i != 0)
                        modInfo += ", ";
                    modInfo += varmods[i].position;
                }
                modInfo += "]</td>";
                modInfo += "</tr>";
	    }
            modInfo += "</table>";
	    modInfo += '</div>';
	}

	$(getElementSelector(container, elementIds.modInfo)).append(modInfo);
    }

    //---------------------------------------------------------
    // ANTIC INFO
    //---------------------------------------------------------
    function showAnticInfo (container) {

        var antic = container.data("anticPrecursor");
        if(labelImmoniumIons(container))
        {
            antic += container.data("anticImmonium");
        }
        if(labelReporterIons(container))
        {
            antic += container.data("anticReporter");
        }
        // selected ions
        var selectedIons = getSelectedIonTypes(container);
        var ionSeriesAntic = container.data("ionSeriesAntic");
        for(var i = 0; i < selectedIons.length; i += 1)
        {
            var sion = selectedIons[i];
            antic += ionSeriesAntic[sion.type][sion.charge];
        }

        // user reporter ions
        if(container.data("anticUserReporterIons")) {
            antic += container.data("anticUserReporterIons");
        }

	var percent = 100*round(antic / container.data("totalIonCurrent"));

	var anticInfo = '<div id="'+getElementId(container, elementIds.anticInfo)+'" style="margin-top:5px;">';
	anticInfo += 'Annotated Ion Current: ';
	anticInfo += "<div><b>"+round(antic)+"</b> ("+percent+"%)</div>";
	anticInfo += '</div>';

	$(getElementSelector(container, elementIds.anticInfo)).replaceWith(anticInfo);
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
	myContent += '<label><input id="'+getElementId(container, elementIds.enableTooltip)+'" type="checkbox">Enable tooltip </label>';
	myContent += '</nobr> ';

        // mass error plot option
        myContent += '<nobr>';
        myContent += '<label><input id="'+getElementId(container, elementIds.massErrorPlot_option)+'" type="checkbox" ';
        if(options.showMassErrorPlot === true)
        {
            myContent += ' checked="checked"';
        }
        myContent += '>Plot mass error </label>';
        myContent += '</nobr>';

	myContent += '<br>';

	$(getElementSelector(container, elementIds.viewOptionsDiv)).append(myContent);
	if(!options.showViewingOptions) {
            $(getElementSelector(container, elementIds.viewOptionsDiv)).hide();
        }
    }


    //---------------------------------------------------------
    // OPTIONS TABLE
    //---------------------------------------------------------
    function makeOptionsTable(container, defaultChargeStates, defaultSelectedIons) {

        var options = container.data("options");

	var myTable = '';
	myTable += '<table cellpadding="2" cellspacing="2"> ';
	myTable += '<tbody> ';

	// Ions
	myTable += '<tr><td class="optionCell"> ';

	myTable += '<b>Ions:</b> ';
	myTable += '<div id="'+getElementId(container, elementIds.ion_choice)+'" style="margin-bottom: 10px"> ';
	myTable += '<!-- a ions --> ';
        myTable += addIonToOptionsTable('a', defaultChargeStates, defaultSelectedIons['a']);
	myTable += '<br/> ';
	myTable += '<!-- b ions --> ';
        myTable += addIonToOptionsTable('b', defaultChargeStates, defaultSelectedIons['b']);
	myTable += '<br/> ';
	myTable += '<!-- c ions --> ';
        myTable += addIonToOptionsTable('c', defaultChargeStates, defaultSelectedIons['c']);
	myTable += '<br/> ';
	myTable += '<!-- x ions --> ';
        myTable += addIonToOptionsTable('x', defaultChargeStates, defaultSelectedIons['x']);
	myTable += '<br/> ';
	myTable += '<!-- y ions --> ';
        myTable += addIonToOptionsTable('y', defaultChargeStates, defaultSelectedIons['y']);
	myTable += '<br/> ';
	myTable += '<!-- z ions --> ';
        myTable += addIonToOptionsTable('z', defaultChargeStates, defaultSelectedIons['z']);
	myTable += '<br/> ';
	myTable += '<span id="'+getElementId(container, elementIds.deselectIonsLink)+'" style="font-size:8pt;text-decoration: underline; color:sienna;cursor:pointer;">[Deselect All]</span> ';
	myTable += '</div> ';

	myTable += '<span style="font-weight: bold;">Neutral Loss:</span> ';
	myTable += '<div id="'+getElementId(container, elementIds.nl_choice)+'"> ';
        var peptide = container.data("options").peptide;
        var idx = 0;
        for (lossKey in peptide.lorikeetPotentialLosses)
        {
            var loss = peptide.lorikeetPotentialLosses[lossKey];
            if(!loss)
                continue;
            if(idx++ != 0)myTable += '<br> ';
            myTable += '<nobr><label> <input type="checkbox" value="'+loss.label()+'" id="'+loss.label()+'"/> ';
            myTable += loss.htmlLabel();
            myTable += '</label></nobr> ';
        }
        for(var lossKey in peptide.customPotentialLosses)
        {
            var loss = peptide.customPotentialLosses[lossKey];
            if(!loss)
                continue;
            if(idx++ != 0)myTable += '<br> ';
            myTable += '<nobr><label> <input type="checkbox" value="'+loss.label()+'" id="'+loss.label()+ '" checked = "checked"/> ';
            myTable += loss.htmlLabel();
            myTable += '</label></nobr> ';
        }
	myTable += '</div> ';

        // Immonium ions
        myTable+= '<label><input type="checkbox" value="true" ';
        if(options.labelImmoniumIons == true)
        {
            myTable+=' checked="checked"';
        }
        myTable+= ' id="'+getElementId(container, elementIds.immoniumIons)+'"/><span style="font-weight:bold;">Immonium ions</span></label>';

        // Reporter ions
        myTable += "<br/>";
        myTable+= '<label><input type="checkbox" value="true" ';
        if(options.labelReporters == true)
        {
            myTable+=' checked="checked"';
        }
        myTable+= ' id="'+getElementId(container, elementIds.reporterIons)+'"/><span style="font-weight:bold;">Reporter ions</span></label>';

        // Unfragmented precursor ions
        myTable += "<br/>";
	myTable+= '<label><input type="checkbox" value="true" ';
        if(options.labelPrecursorPeak === true)
	{
	    myTable+=' checked="checked"';
	}
        myTable+= ' id="'+getElementId(container, elementIds.labelPrecursor)+'"/><span style="font-weight:bold;">Precursor ions</span></label>';

	myTable += '</td> </tr> ';

	// mass type
	myTable += '<tr><td class="optionCell"> ';
	myTable += '<div> Frag. Mass Type:<br/> ';
	myTable += '<nobr> ';
	myTable += '<label><input type="radio" name="'+getRadioName(container, "massTypeOpt")+'" value="mono"';
        if(options.fragmentMassType == 'mono')
            myTable += ' checked = "checked" ';
        myTable += '/><span style="font-weight: bold;">Mono</span></label> ';
	myTable += '<label><input type="radio" name="'+getRadioName(container, "massTypeOpt")+'" value="avg"';
        if(options.fragmentMassType == 'avg')
            myTable += ' checked = "checked" ';
        myTable += '/><span style="font-weight: bold;">Avg</span></label> ';
	myTable += '</nobr> ';
	myTable += '</div> ';
        // mass tolerance
	myTable += '<div style="margin-top:10px;"> ';
	myTable += '<nobr>';
        myTable += 'Mass Tol: <input id="'+getElementId(container, elementIds.massError)+'" type="text" value="'+options.massError+'" style="width:3em;"/>';
	myTable += '</nobr>';
        myTable += '<br>';
        myTable += '<label><input type="radio" name="'+getRadioName(container, "massErrorUnit")+'" value="' + massErrorTypeTh + '"';
        if(options.massErrorUnit === massErrorTypeTh)
        {
            myTable += ' checked = "checked" ';
        }
        myTable += '/><span style="font-weight: bold;">' + massErrorTypeTh + '</span></label> ';
        myTable += '<label><input type="radio" name="'+getRadioName(container, "massErrorUnit")+'" value="' + massErrorTypePpm + '"';
        if(options.massErrorUnit === massErrorTypePpm)
        {
            myTable += ' checked = "checked" ';
        }
        myTable += '/><span style="font-weight: bold;">' + massErrorTypePpm + '</span></label> ';
        myTable += '</div> ';
	myTable += '<div style="margin-top:10px;" align="center"> ';
	myTable += '<input id="'+getElementId(container, elementIds.update)+'" type="button" value="Update"/> ';
	myTable += '</div> ';
	myTable += '</td> </tr> ';

	// peak assignment method
	myTable += '<tr><td class="optionCell"> ';
	myTable+= '<div> Peak Assignment:<br/> ';
	myTable+= '<label><input type="radio" name="'+getRadioName(container, "peakAssignOpt")+'" value="intense" checked="checked"/><span style="font-weight: bold;">Most Intense</span></label><br/> ';
	myTable+= '<label><input type="radio" name="'+getRadioName(container, "peakAssignOpt")+'" value="close"/><span style="font-weight: bold;">Nearest Match</span></label><br/> ';
        myTable+= '<label><input type="checkbox" value="true" ';
        if(options.peakDetect == true)
        {
            myTable+=' checked="checked"';
        }
        myTable+= ' id="'+getElementId(container, elementIds.peakDetect)+'"/><span style="font-weight:bold;">Peak Detect</span></label>';
	myTable+= '</div> ';
	myTable += '</td> </tr> ';

	// peak labels
	myTable += '<tr><td class="optionCell"> ';
	myTable+= '<div> Peak Labels:<br/> ';
	myTable+= '<label><input type="radio" name="'+getRadioName(container, "peakLabelOpt")+'" value="ion" checked="checked"/><span style="font-weight: bold;">Ion</span></label>';
	myTable+= '<label><input type="radio" name="'+getRadioName(container, "peakLabelOpt")+'" value="mz"/><span style="font-weight: bold;">m/z</span></label><br/>';
	myTable+= '<label><input type="radio" name="'+getRadioName(container, "peakLabelOpt")+'" value="none"/><span style="font-weight: bold;">None</span></label> ';
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
        if(!options.showOptionsTable) {
            $(getElementSelector(container, elementIds.optionsTable)).hide();
        }
    }

    function addIonToOptionsTable(ionLabel, charges, selected)
    {
        if(!selected) selected = [];
        var ionRow = "";
        ionRow += '<nobr> ';
        ionRow += '<span style="font-weight: bold;">' + ionLabel + '</span> ';
        for (var i = 0; i < charges.length; i += 1)
        {
            var id = ionLabel + "_" + charges[i];
            var checked = (selected[i] && selected[i] == 1) ? ' checked="checked"' : "";
            ionRow += '<input type="checkbox" value="' + charges[i] + '" id="' + id + '" ' + checked + '/>' + charges[i] + '<sup>+</sup> ';
        }
        ionRow += '</nobr> ';
        return ionRow;
    }


})(jQuery);
