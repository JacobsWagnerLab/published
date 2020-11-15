%{
-------------
  PURPOSE:   
-------------
The basic purpose of this class is to create a pipeline to go from a map,
or network to the display of enrichment results over this map/network
with a minimal amount of intermediate steps, and in the most robust
possible way (that I could think of).

Define the mapExplorer class to compute local enrichments over a network
The nodes in the network are characterized by a set of features that are
used to generate the network and each node is characterized by a set of
attribute values for which we aim to calculate enrichement probabilities

The idea of the code initially tsarted to take shape with Bruno Beltran
when we started to use tSNE algortihm. However, the major part of the
algorithm is based on the publication by Anastasia Baryshnikova
      (Systematic %functional annotation and visualization of biological
       networks, (2016) Cell systems 2 (6), 412-421)
describing a very similar idea for networks in general. I
implemented a modified version of this algorithm in this MATLAB class 
mapExplorer to couple it to tSNE maps and all sorts of features and 
attributes.

JUSTIFY ENRICHMENT FOR FEATURES

-------------
 DEPENDENCES
-------------
*kde2dv2.m
      available at:
      http://www.mathworks.com/matlabcentral/fileexchange/
          17204-kernel-density-estimation


*compute_mapping.m from the  Matlab Toolbox for Dimensionality Reduction
      reference for tSNE & website for the Matlab toolbox:
      L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional
      Data Using t-SNE. Journal of Machine Learning Research
      9(Nov):2579-2605, 2008.

      https://lvdmaaten.github.io/drtoolbox/
      
*fdrbh.m


INITIALIZATION
MmapObj = mapExplorer(map, features, featureNames, attributes,...
                attributesNames, attributesStyle, distPct)

Explain


FUNCTION CALLS
__________________
mapExplorer

Syntax:
MmapObj = mapExplorer(map, features, featureNames, attributes,...
                attributesNames, attributesStyle, distPct)

Inputs: 
  map               n x 2 array defining teh (x,y) positions of the
                    nodes in the map/network. If empty, a tSNE map 
                    based on the feature input will be calculated.
  features          n x m double array. Required input.
  featureNames      m x 1 or 1 x m cell array of strings giving the
                    names of the m different features.
                    If empty, a default cell array of 'FeatureXXX' will
                    be created
  attributes        n x k double array. Required input
  attributesNames   k x 1 or 1 x k cell array of strings giving the
                    names of the k different attributes.
                    If empty, a default cell array of 'AttributeXXXX' 
                    will be created
  attributesStyle   Type of attributes used - define how enrichments
                    are evaluated
                    - 'binary'
                    - 'continuous'
  distPct           Define the percentile of the pairwise node
                    distances distribution to be used as a threshold.
                    Typically between 1 and 5%. See Baryshnikova (2016)

Output:
  MmapObj         An object with all the fields mentionned above filled

=========================

  PUBLIC CLASS METHODS

=========================

setAttributes

Syntax:
mapObj = setAttributes(mapObj, attributes, attributesNames)

Inputs: 
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  attributes      n x k double array of n instances of measurements for k
                  attributes.
  attributesNames k x 1 (or 1 x k) cell array of strings giving the
                  names of the k different attributes.

Output:
  mapObj         An object with all the fields mentionned above filled
                 (and specifically the attributes and attributesNames
                 fields)

_________________________
setFeatures

Syntax:
mapObj = setFeatures(mapObj, features, map, featureNames)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  features        n x m double array of n instances of measurements for k
                  features.
  featureNames    m x 1 (or 1 x m) cell array of strings giving the
                  names of the m different features.


Output
  mapObj         An object with all the fields mentionned above filled
                 (and specifically the features and featureNames
                 fields)

_________________________
setAttributesNames

Syntax:
mapObj = setAttributesNames(mapObj, attributesNames)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  attributesNames k x 1 (or 1 x k) cell array of strings giving the
                  names of the k different attributes.

Output
  mapObj         An object with all the fields mentionned above filled
                 (and specifically the  attributesNames fields)

_________________________
setFeaturesNames

Syntax:
mapObj = setFeaturesNames(mapObj, featuresNames)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  featureNames    m x 1 (or 1 x m) cell array of strings giving the
                  names of the m different features.

Output
  mapObj         An object with all the fields mentionned above filled
                 (and specifically the  featureNames fields)

_________________________
setDistanceThreshold

Syntax:
mapObj = setDistanceThreshold(mapObj, distance_value)

Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  distance_value  Value from the ]0, 1] interval corresponding to the
                  percentile of the pairwise inter-node distances to
                  calculate the distance threshold defining the radius of
                  the neighborhood in which enrichments are calculated.
                  0.05 means the 5th percentile of the distance
                  distribution.


Output
  mapObj         An object with all the fields mentionned above filled


_________________________
plotDensityEnvelope

Syntax:
hf = plotDensityEnvelope(mapObj)


Input:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).

Output
  hf              handle to the figure displaying the density contour of
                  the tSNE map/network

_________________________
getDistances

Syntax:
dst = getDistances(mapObj)


Input:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).

Output
  dst             percentile of the pairwise inter-node distances defined
                  by the input variable 'distPct', or via the method
                  'setDistanceThreshold'

_________________________
getNeighbors

Syntax:
neigh = getNeighbors(mapObj)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).

Output
  neigh           n x n logical array specifying which vector (line in
                  the feature array) is neighbor, or in other word within
                  a given neighborhood, with which other feature vector
                  in the tSNE map/network.
                  0/false:	not neighbors
                  1/true:     neighbors

_________________________
getRadius

Syntax:
r = getRadius(mapObj)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).

Output
  r               distance value corresponding to the percentile of the
                  pairwise inter-node distances defined by the input
                  variable 'distPct', or via the method
                  'setDistanceThreshold'

_________________________
getDensity

Syntax:
density = getDensity(mapObj)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).

Output
  density         256 x 256 double array defining the probability density
                  of the nodes in teh tSNE map/network. The density
                  probability is derived via kernel density estimator
                  function. See kde2d.m as a dependent function.

_________________________
getMapEnvelope

Syntax:
envelope = getMapEnvelope(mapObj)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).

Output
  envelope        (x,y) two columns array defining the contour of the
                  probability density matrix at the level specified in
                  the function (typically around 1e-4 for ~1,000 points)

_________________________
getPval

Syntax:
pval = getPval(mapObj, flag)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  flag            Either 'features' or attributes'. Defines which of the
                  two data sets are used to calculate local enrichments.

Output
  pval            n x m double array of the p-values for the
                  hypergeometric test for the neighborhood around each
                  node in the tSNE map/network and for each feature or
                  attribute (depending on the value of the input vairable
                  'flag').

_________________________
getQval

Syntax:
qval = getQval(mapObj, flag)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  flag            Either 'features' or attributes'. Defines which of the
                  two data sets are used to calculate local enrichments.

Output
  qval            n x m double array of the q-values for the
                  hypergeometric test for the neighborhood around each
                  node in the tSNE map/network and for each feature or
                  attribute (depending on the value of the input vairable
                  'flag'). The p-values are corrected based on the false
                  discovery rate as defined by Benjamini & Yekutieli
                  (2001) procedure. The procedure is encoded in the
                  fdrbh.m

_________________________
plotEnrichedNeighborhoods

Syntax:
hfig = plotEnrichedNeighborhoods(mapObj, flag, qthresh)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  flag            Either 'features' or attributes'. Defines which of the
                  two data sets are used to calculate local enrichments.
  qthresh         FDR (or q-value) threshold signifying the significance
                  threshold to consider a neighborhood to be
                  statistically enriched (typically 0.05 or 0.01).

Output
  hfig            Cell array of handles to the figures displaying the
                  density contour of the tSNE map/network and
                  highlighting the location of the enriched neighborhoods
                  on this map.
                  One figure for each feature or attribute is generated
                  in the maximum of 25 figures. For each feature or
                  attribute (depending on the value of the input vairable
                  'flag') all the enriched neiborhood are highlighted by
                  a point with a color intensity incressing with the
                  -log(q-value) value scaled to the interval
                  [-log10(qtresh) , -log10(10^-8)]

_________________________
plotFeatureMap

Syntax:
hfig = plotFeatureMap(mapObj, threshScore)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  threshScore     Minimal absolute score used to select the nodes to plot
                  within the density probability contour, i.e. setting
                  threshSore to 3 will allow the display of all the nodes
                  in the map associated with a feature score above 3 or
                  below -3.

Output
  hfig            Cell array of handles to the figures displaying the
                  density contour of the tSNE map/network and
                  highlighting the location of the nodes associated with
                  a feature score below or above threshScore.
                  One figure per feature is generated.

_________________________
plotAttributeMap

Syntax:
hfig = plotAttributeMap(mapObj, qthresh)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  qthresh         FDR (or q-value) threshold below which an attribute is
                  deemed to be significantly enriched (typically 0.05 or
                  0.01).

Output
  hfig            Handle to the figures displaying the
                  density contour of the tSNE map/network and
                  highlighting the location of the enriched neighborhoods
                  One figure is generated for all attribute. The density
                  of the purple color relates to the value of the q-value
                  scaled to the [-log10(qtresh) , -log10(10^-8)] interval

_________________________
plotAttributeHeatmap

Syntax:
[hfig, H, Hgo] = plotAttributeHeatmap(mapObj, qthresh)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).
  qthresh         FDR (or q-value) threshold below which an attribute is
                  deemed to be significantly enriched (typically 0.05 or
                  0.01).

Output
  hfig            Handle to the figures displaying the clustergram of the
                  average "phenotypes", or average feature vector 
                  (averaged over the nodes present in the enriched
                  neighborhoods for a given attribute).
                  The clustergram is generated with the dependent
                  function clusterMap.m
                  The color scale is set to the [-6 6] range by default.

_________________________
getContours

Syntax:
contours = getContours(mapObj)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).

Output
  contours        cell array of vectors defining the (x,y) positions for
                  all the enriched areas separated by feature and either
                  under- and over-represented areas of teh tSNE
                  map/network


=========================

 PRIVATE CLASS METHODS

=========================

_________________________
featureEnrichement

Syntax:
mapObj = featureEnrichement(mapObj)

Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).

Output
  mapObj          mapExplorer object filled with the field pvalFeatures.
                  n x m array of p-values calculated for the enrichemnt
                  of each feature around the neighborhoods centered on
                  each node of the tSNE map/network.

_________________________
attributeEnrichement

Syntax:
mapObj = attributeEnrichement(mapObj)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method). 

Output
  mapObj          mapExplorer object filled with the field pvalAttributes.
                  n x k array of p-values calculated for the enrichemnt
                  of each attribute around the neighborhoods centered on
                  each node of the tSNE map/network.

_________________________
correctFDRfeatures

Syntax:
mapObj = correctFDRfeatures(mapObj, varargin)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method). 

Output
  mapObj          mapExplorer object filled with the field qvalFeatures.
                  n x m array of q-values calculated by correcting the
                  pvalues (pvalFeatures) via the control of the false
                  discovery rate as defined by Benjamini & Yekutieli.
                  See the dependent function fdrbh for this correction

_________________________
correctFDRattributes

Syntax:
mapObj = correctFDRattributes(mapObj, varargin)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method). 

Output
  mapObj          mapExplorer object filled with the field qvalAttributes.
                  n x k array of q-values calculated by correcting the
                  pvalues (pvalAttributes) via the control of the false
                  discovery rate as defined by Benjamini & Yekutieli.
                  See the dependent function fdrbh for this correction

_________________________
getFeatureContours

Syntax:
mapObj = getFeatureContours(mapObj, varargin)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).  

Output
  mapObj          mapExplorer object filled with the field contours.
                  2D cell array of 3D arrays  defining the contours
                  of the areas of the tSNE map/network where attributes
                  are enriched.
                  The cell array 2D:

                  The double array 3D:

                  This is an outdated and ugly method.


The last two functions are separated for simplicity of use: A single 
output for each is required, only. The drawback is that you need to call 
both functions if you want both Xbins and Ybins. See below for 
explanations about the use and the rationale of these two functions.

_________________________
getXbins

Syntax:
Xbins = getXbins(mapObj)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).  

Output
  Xbins           Array of the x-axis binning array used to compute the
                  probability density map (256 entry vector)

_________________________
getYbins

Syntax:
Ybins = getYbins(mapObj)


Inputs:
  mapObj          mapExplorer object (at the very least initialized with
                  the mapExplorer method).  

Output
  Ybins           Array of the y-axis binning array used to compute the
                  probability density map (256 entry vector)




===========================
Manuel Campos
Jacobs-Wagner lab
August 2016
%}

classdef mapExplorer < handle
    
    properties
        features
        featureNames
        map
        attributes
        attributesNames
        attributesStyle
        distPct
    end
    
    properties (Access = private)
        mapEnvelope 
        mapDensity 
        mapXbins 
        mapYbins 
        distances 
        radius 
        neighbors 
        pvalFeatures 
        qvalFeatures 
        enrNeighFeatures 
        pvalAttributes 
        qvalAttributes 
        enrNeighAttributes 
        contours 
    end
    
    methods
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %       CLASS CONSTRUCTOR
        %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initialize the mapExplorer object
        % Incorporate the input provided
        % Proceed to calculate a tSNE map (if necessary) based on 'features'
        % Proceed to calculate enrichments of the 'attributes'
        
        function mapObj = mapExplorer(map, features, featureNames, attributes,...
                attributesNames, attributesStyle, distPct)
            % develop robustness,mostly vis a vis inputs
            if nargin<2 || isempty(features) || isempty(attributes)
                disp(['You need to define at the very least both',...
                    ' "features" and "attributes" properties']);
                return
            end
            % Initialize the object and fill it in based on the input
            % variables available.
            mapObj.map = map;
            mapObj.features = features;
            mapObj.attributes = attributes;
            % If the provided map is empty ([]), calculate a 2D tSNE map
            % from the n x m features
            if isempty(map)
                disp('Computing a 2D tSNE map based on the features');
                mapObj.map = compute_mapping(features, 'tSNE', 2);
            end
            % Create default attribute names cell array of string with a 4
            % digits number. If you have more than 9,999 features, you need
            % (a) to provide attribute names
            % (b) change the definition of teh automatic numbering to %1.5u
            % (c) do not use this class as it is not optimized for such
            %  large datasets and will likely use up all your RAM :D
            if isempty(attributesNames)
                mapObj.attributesNames = cell(size(attributes,2),1);
                for ii = 1:size(attributes,2)
                    mapObj.attributesNames{ii} = ['Attribute',num2str(ii,'%4.4u')];
                end
            else
                mapObj.attributesNames = attributesNames;
            end
            % Similarly, create default attribute names cell array of 
            % string with a 3 digits number. For more than 999 features,
            % see comment above.
            if isempty(featureNames)
                mapObj.featureNames = cell(size(features,2),1);
                for ii = 1:size(features,2)
                    mapObj.featureNames{ii} = ['Feature',num2str(ii,'%3.3u')];
                end
            else
                mapObj.featureNames = featureNames;
            end
            % If the input variable 'attributesStyle' is left empty ([]),
            % try to guess waht would be the best attribute style depending
            % on the attribute data themselves.
            % If there are only two unique values in the attribute array,
            % associate the attributeStyle variable with 'binary'.
            % If there are more than 2 unique values, make it 'continuous'
            %
            % ALSO: default distance threshold defining the siz eof the
            % neighborhood ahere enrichments are calculated is set to 5% by
            % default
            if nargin<6
                mapObj.attributesStyle = 'binary';
                mapObj.distPct = 0.5;
                if length(unique(attributes(:)))>2
                    mapObj.attributesStyle = 'continuous';
                end
            elseif isempty(attributesStyle)
                mapObj.attributesStyle = 'binary';
                if length(unique(attributes(:)))>2
                    mapObj.attributesStyle = 'continuous';
                end
            else
                mapObj.attributesStyle = attributesStyle;
            end
            if nargin==7
                mapObj.distPct = distPct;
            end
            
            % Check inputs coherence
            if ~isequal(size(mapObj.features,2),length(mapObj.featureNames))
                disp(['The number of lines in the properties "features"',...
                    ' and "featuresNames" must match']);
                return
            end
            if ~isequal(size(mapObj.attributes,2),length(mapObj.attributesNames))
                disp(['The number of lines in the properties "attributes"'...
                    ' and "attributesNames" must match']);
                return
            end
            if ~isempty(mapObj.map) && (~isequal(size(mapObj.map,1),size(mapObj.features,1)) || ...
                    ~isequal(size(mapObj.map,1),size(mapObj.attributes,1)))
                disp(['The number of lines in the properties "map",'...
                    ' "features" and "attributes" must match']);
                return
            end
            if ~isempty(mapObj.map) && (~isequal(size(mapObj.map,2),2))
                disp('The similarity map should be 2-dimensional.');
                fprintf(2,...
                    'error in constructor: property size size(mapObj.map,2) ~= 2');
                return
            end
            if ~isempty(mapObj.distPct) && (mapObj.distPct<=0 || mapObj.distPct>=100)
                disp(['The property "distPct" value is a percentile',...
                    ' of the pairwise distances distributions.']);
                disp('The input value should be in the range [0 1]');
                return
            end
            % fill all possible properties here
            mapObj.distances = squareform(getDistances(mapObj));
            mapObj.radius = getRadius(mapObj);
            mapObj.neighbors = getNeighbors(mapObj);
            mapObj = featureEnrichement(mapObj);
            mapObj = attributeEnrichement(mapObj);
            mapObj = correctFDRfeatures(mapObj, 0.05);
            mapObj = correctFDRattributes(mapObj, 0.05);
            mapObj = getFeatureContours(mapObj, 0.05);
            plotEnrichedNeighborhoods(mapObj, 'attributes', 0.05);
            %Add a line to save figures?
        end
    
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %       PUBLIC CLASS METHODS
        %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function mapObj = setAttributes(mapObj, attributes, attributesNames)
%             if any(isnan(attributes(:)))
%                 disp('NaN values not supported');
%                 return
%             else
                mapObj.attributes = attributes;
                mapObj = attributeEnrichement(mapObj);
                mapObj = correctFDRattributes(mapObj, 0.05);
%             end
            if nargin<3 || isempty(attributesNames)
                mapObj.attributesNames = cell(size(attributes,2),1);
                for ii = 1:size(attributes,2)
                    mapObj.attributesNames{ii} = ['Attribute',num2str(ii,'%4.4u')];
                end
            elseif size(attributesNames,1) == size(attributes,2)
                mapObj.attributesNames = attributesNames;
            else
                fprintf(2,'The number of lines in the properties "attributes" and "attributesNames" must match\n');
                return
            end
            if ~isequal(size(mapObj.features,1),size(mapObj.attributes,1))
                fprintf('\n');
                fprintf(2,'The number of lines in the properties "attributes" and "features" must match\n');
                fprintf(2,'An update of the "features" property is necessary\n');
            end
        end
        
        function mapObj = setFeatures(mapObj, features, map, featureNames)
            if any(isnan(features(:)))
                disp('NaN values not supported');
                return
            else
                if nargin<3 || isempty(map)
                    map = compute_mapping(features, 'tSNE', 2);
                end
                mapObj.map = map;
                mapObj.features = features;
                mapObj = featureEnrichement(mapObj);
                mapObj = correctFDRfeatures(mapObj, 0.05);
                mapObj = getFeatureContours(mapObj, 0.05);
            end
            if nargin<4 || isempty(featureNames)
                mapObj.featureNames = cell(size(features,2),1);
                for ii = 1:size(features,2)
                    mapObj.featureNames{ii} = ['Feature',num2str(ii,'%3.3u')];
                end
            elseif length(featureNames) == size(features,2)
                mapObj.featureNames = featureNames;
            else
                disp('The number of lines in the properties "features" and "featuresNames" must match');
                return
            end
            if ~isequal(size(mapObj.features,1),size(mapObj.attributes,1))
                fprintf('\n');
                fprintf(2,'The number of lines in the properties "attributes" and "features" must match\n');
                fprintf(2,'An update of the "attributes" property is necessary\n');
            end
        end
        
        function mapObj = setAttributesNames(mapObj, attributesNames)
            if iscell(attributesNames) && ...
                    ischar(attributesNames{1}) &&...
                    size(attributesNames,1)==size(mapObj.attributes,2)
                mapObj.attributesNames = attributesNames;
            else
                disp('Attributes names must be set as a cell array of characters');
                disp('with the same number of lines as the property "attributes"');
            end
        end
        
        function mapObj = setFeaturesNames(mapObj, featuresNames)
            if iscell(featuresNames) && ...
                    ischar(featuresNames{1}) &&...
                    length(featuresNames)==size(mapObj.features,2)
                mapObj.featureNames = featuresNames;
            else
                disp('Features names must be set as a cell array of characters');
                disp('with the same number of columns as the property "features"');
            end
        end
        
        function mapObj = setDistanceThreshold(mapObj, value)
            if value>0 && value<100
                mapObj.distPct = value;
                mapObj.distances = squareform(getDistances(mapObj));
                mapObj.radius = getRadius(mapObj);
                mapObj.neighbors = getNeighbors(mapObj);
                mapObj = featureEnrichement(mapObj);
                mapObj = attributeEnrichement(mapObj);
                mapObj = correctFDRfeatures(mapObj, 0.05);
                mapObj = correctFDRattributes(mapObj, 0.05);
                mapObj = getFeatureContours(mapObj, 0.05);
            else
                disp('The input value is a percentile of the pairwise distances distributions.');
                disp('The input value should be in the range [0 100]');
            end
        end
        
        function  hf = plotDensityEnvelope(mapObj) 
            [~,density,X,Y]=kde2dv2(mapObj.map,2^8);
            bnd = boundary(mapObj.map(:,1),mapObj.map(:,2),0.6);
            IN = inpolygon(X(:),Y(:),mapObj.map(bnd,1),mapObj.map(bnd,2));
            msk = density;
            msk(IN) = 1;
            msk(~IN) = 0;
            se = strel('disk',3);
            mskdil = imdilate(msk,se);
            figure;envelope = contour(X,Y,mskdil,[0.1 0.1]);close(gcf);

            hf = figure;
            pcolor(X,Y,density);shading flat;colormap(flipud(bone));
            hold on;plot(envelope(1,2:end),envelope(2,2:end),...
                '-','linewidth',1.5);
            set(gca,'xlim',[min(mapObj.map(:,1))-10, max(mapObj.map(:,1))+10],...
                'ylim',[min(mapObj.map(:,2))-10, max(mapObj.map(:,2))+10]);
            axis equal;axis off;
        end
        
        function dst = getDistances(mapObj)
            dst = pdist(mapObj.map,'euclidean');
        end
        
        function r = getRadius(mapObj)
            if isempty(mapObj.distPct)
                mapObj.distPct = 0.5;
                fprintf('\n');
                disp('Calculating the radius based on the default (5%) percentile of pairwise distances');
                disp('Use "setDistanceThreshold" method to modulate this value');
            end
            dst = getDistances(mapObj);
            r = prctile(dst, mapObj.distPct); 
        end
        
        function neigh = getNeighbors(mapObj)
            dst = getDistances(mapObj);
            r = getRadius(mapObj);
            dst = squareform(dst);
            neigh = dst <= r;
        end
        
        function density = getDensity(mapObj)
            [~,density,~,~]=kde2dv2(mapObj.map,2^8);
        end
        
        function envelope = getMapEnvelope(mapObj)
            [~,density,X,Y]=kde2dv2(mapObj.map,2^8);
            bnd = boundary(mapObj.map(:,1),mapObj.map(:,2),0.6);
            IN = inpolygon(X(:),Y(:),mapObj.map(bnd,1),mapObj.map(bnd,2));
            msk = density;
            msk(IN) = 1;
            msk(~IN) = 0;
            se = strel('disk',3);
            mskdil = imdilate(msk,se);
            figure;
            envelope = contour(X,Y,mskdil,[0.1 0.1]);close(gcf);
%             envelope = contour(X,Y,density,[5e-5 5e-5],'-k');close(gcf);
        end

        function pval = getPval(mapObj, flag)
            if strcmpi('features', flag)
                pval = mapObj.pvalFeatures;
            elseif strcmpi('attributes', flag)
                pval = mapObj.pvalAttributes;
            else
                disp('Error in getPval method.');
                disp('Second input: {"features" | "attributes"}');
            end
        end
        
        function qval = getQval(mapObj, flag)
            if strcmpi('features', flag)
                qval = mapObj.qvalFeatures;
            elseif strcmpi('attributes', flag)
                qval = mapObj.qvalAttributes;
            else
                disp('Error in getQval method.');
                disp('Second input: {"features" | "attributes"}');
            end
        end
        
        function hfig = plotEnrichedNeighborhoods(mapObj, flag, qthresh)
            % Selected feature of interest based on FDR corrected p-values
            if nargin==1
                flag = 'attributes';
                qthresh = 0.05;
                disp('Plotting the neighborhoods enriched in some attributes (default)');
                disp('Using the default FDR level (5%)');
            end
            if nargin>1 && ~any(strcmpi(flag,{'features';'attributes'}))
                flag = 'attributes';
                disp('Plotting the neighborhoods enriched in attributes (default)');
            end
            if nargin<3 || isempty(qthresh)
                qthresh = 0.05;
            elseif qthresh<=0 || qthresh>=1
                qthresh = 0.05;
                fprintf('\n');
                fprintf(2,'The q-value threshold must be in the range [0 1]\n');
                fprintf('Using the default FDR level (5%%)\n');
            else
                fprintf('\n');
                fprintf('Display enriched %s at an FDR level of %f\n', flag, qthresh);
            end
            % Prepare radius legend
            r = getRadius(mapObj);
            th = linspace(0,2*pi,100);
            xc = r.*cos(th) + max(mapObj.map(:,1)) - r/2;
            yc = r.*sin(th) + min(mapObj.map(:,2)) + r/2;
            
            names = mapObj.attributesNames;
            dataType = mapObj.attributesStyle;
            if strcmpi(flag,'features');
                names = mapObj.featureNames ;
                dataType = 'continuous';
            end

            qval = getQval(mapObj, flag);
            enrNeigh = qval <= qthresh;
            nn = find(sum(enrNeigh(:,:,1))>0);
            if strcmpi(dataType,'continuous')
                nn = unique([find(sum(enrNeigh(:,:,1))>0), find(sum(enrNeigh(:,:,2))>0)]);
            end
            hfig = {};
            if ~isempty(nn)
            hfig = cell(length(nn),1);
            for ii = nn(1:min([length(nn),25])) % Plot top 25 features max
                % Scale the qvalues
                S_high = -log10(qval(:,ii,1) + eps)./8;%(-log10(min(enrMap.qval(:))+eps));
                S_high(S_high>1) = 1; % Saturate scores at q-values = 10^-8 or below
                S_low  = -log10(qval(:,ii,2) + eps)./8;%(-log10(min(enrMap.qval(:))+eps));
                S_low(S_low>1) = 1; % Saturate scores at q-values = 10^-6 or below

                c = ones(length(mapObj.map),3);
                % Enriched areas locations and colors
                ix_high = qval(:,ii,1)<qthresh;
                chigh = repmat([0 .8 1],sum(ix_high),1);
                chigh = 1 - bsxfun(@times,chigh,S_high(ix_high));
                c(ix_high,:) = chigh;

                if strcmpi('continuous',dataType)
                    ix_low = qval(:,ii,2)<qthresh;
                    clow = repmat([1 .6 0],sum(ix_low),1);
                    clow = 1 - bsxfun(@times,clow,S_low(ix_low));
                    c(ix_low,:) = clow;
                end

                % Sort data by color intensity
                [csum,ix] = sort(sum(c,2),'descend');
                x = mapObj.map(ix,1);
                y = mapObj.map(ix,2);
                c = c(ix,:);
                hfig{ii} = figure; hold on; %set(gcf,'Color','k');
                scatter(x(csum<2.99), y(csum<2.99), 120, c(csum<2.9,:), 'filled');
%                 envelope = getMapEnvelope(mapObj);
%                 plot(envelope(1,2:end),envelope(2,2:end),'-',...
%                     'linewidth',1.5);
                [~,density,X,Y]=kde2dv2(mapObj.map,2^8);
                contour(X,Y,density,[1e-5 1e-5],'-k','linewidth',0.8);hold on;
                set(gca,'xlim',[min(x)-10, max(x)+10],'ylim',[min(y)-10, max(y)+10]);
                title(names{ii});
                plot(xc,yc,':k');
                axis equal;axis off;
                savePath = [];%'Y:\Manuel Campos\Publications\Keio\Figures\tmp\xMor20160831';%
                if ~isempty(savePath)
                print(hfig{ii},'-depsc','-r300',...
                    [savePath,'\xM_F',num2str(ii,'%2.2u')]); % 'Y:\Manuel Campos\Publications\Keio\Figures\fig\Fig3
                end
            end
            else
                disp('No enrichment detected');
            end
        end
        
        function hfig = plotFeatureMap(mapObj, threshScore)
            if nargin == 1
                threshScore = 3;
            end
            scores = mapObj.features;
            scores(abs(scores) < threshScore) = 0;
            scores = scores./10;
            scores = min(max(scores,-1),1);
            map = mapObj.map;
            [~,density,X,Y]=kde2dv2(mapObj.map,2^8);
            hfig = cell(size(mapObj.features,2),1);
            for ii = 1:length(hfig)
                hfig{ii} = figure;hold on;
                S = scores(:,ii);
                
                ix_pos = S>=0;
                cpos = repmat([0 .8 1],sum(ix_pos),1);
                cpos = 1 - bsxfun(@times,cpos,S(ix_pos));
                c(ix_pos,:) = cpos;

                ix_neg = S<0;
                cneg = repmat([1 .6 0],sum(ix_neg),1);
                cneg = 1 - bsxfun(@times,cneg,-S(ix_neg));
                c(ix_neg,:) = cneg;

                [csum,ix] = sort(sum(c,2),'descend');
                x = map(ix,1);
                y = map(ix,2);
                c = c(ix,:);
                
                contour(X,Y,density,[3e-5 3e-5],'-k','linewidth',0.8);hold on;
                scatter(x(csum<2.99), y(csum<2.99), 30, c(csum<2.99,:), 'filled');

                set(gca,'xlim',[min(map(:,1))-10, max(map(:,1))+10],...
                    'ylim',[min(map(:,2))-10, max(map(:,2))+10]);
                axis equal; axis off
                title({mapObj.featureNames{ii};...
                    ['|threshold| = ',num2str(threshScore)]});
            end
        end
        
        function hfig = plotAttributeMap(mapObj, qthresh)
            % Selected feature of interest based on FDR corrected p-values
            if nargin==1 || isempty(qthresh)
                qthresh = 0.05;
                disp('Using the default FDR level (5%)');
            elseif qthresh<=0 || qthresh>=1
                qthresh = 0.05;
                fprintf('\n');
                fprintf(2,'The q-value threshold must be in the range [0 1]\n');
                fprintf('Using the default FDR level (5%%)\n');
            else
                fprintf('\n');
                disp(['Display enriched features at an FDR level of ',num2str(qthresh)]);
            end
            % Prepare radius legend
            r = getRadius(mapObj);
            th = linspace(0,2*pi,100);
            xc = r.*cos(th) + max(mapObj.map(:,1)) - r/2;
            yc = r.*sin(th) + min(mapObj.map(:,2)) + r/2;
            
            dataType = 1;
            if strcmpi('continuous',mapObj.attributesStyle)
                dataType =2;
            end
            qval = getQval(mapObj, 'attributes');
            enrNeigh = qval <= qthresh;
            for jj = 1:dataType
                nn = find(sum(enrNeigh(:,:,jj))>0);
                if ~isempty(nn)
                    clmap = repmat([.65 0.1 .75],sum(sum(enrNeigh(:,:,jj))>0),1);
                    c = ones(length(mapObj.map),3);
                    kk = 1;
                    for ii = nn(1:min([length(nn),50])) % Plot top 20 features max
                        % Scale the qvalues
                        Scores = -log10(qval(:,ii,jj) + eps)./8;%(-log10(min(qval(:))+eps));
                        Scores(Scores>1) = 1; % Saturate scores at q-values = 10^-8 

                        % Enriched areas locations and colors
                        ix_enr = qval(:,ii,jj) < qthresh;
                        c_enr = repmat(1 - clmap(kk,:),sum(ix_enr),1);
                        c_enr = 1 - bsxfun(@times,c_enr,Scores(ix_enr));
                        c(ix_enr,:) = c_enr;
                        kk = kk+1;
                    end

                    % Sort data by color intensity
                    [csum,ix] = sort(sum(c,2),'descend');
                    x = mapObj.map(ix,1);
                    y = mapObj.map(ix,2);
                    c = c(ix,:);
                    hfig = figure;
                    scatter(x(csum<2.99), y(csum<2.99), 120, c(csum<2.99,:), 'filled');hold on;
                    [~,density,X,Y]=kde2dv2(mapObj.map,2^8);
                    contour(X,Y,density,[5e-6 5e-6],'-k','linewidth',0.8);hold on;
                    set(gca,'xlim',[min(mapObj.map(:,1))-10, max(mapObj.map(:,1))+10],...
                        'ylim',[min(mapObj.map(:,2))-10, max(mapObj.map(:,2))+10]);
                    title('Global enrichment map');
                    plot(xc,yc,':k');
                    axis equal;axis off;
                else
                    disp('No enrichment detected');
                end
            end
        end
        
        function [hfig, H, Hgo] = plotAttributeHeatmap(mapObj, qthresh)
            
            dim3 = 1; % 1 - enriched | 2 - empoverished
            neigh = getNeighbors(mapObj);
            qv = getQval(mapObj, 'attributes');
            % For each term with at least one enriched neighborhood,
            % get all neighbors and count them only once
            sigNeigh = find(sum(qv(:,:,dim3)<=0.05)>0);
            GOind = zeros(length(sigNeigh),length(neigh));
            for ii=1:length(sigNeigh)
                locNeigh = neigh(qv(:,sigNeigh(ii),dim3)<=qthresh,:);
                if size(locNeigh,1)==1
                    GOind(ii,:) = locNeigh;
                else
                    GOind(ii,:) = sum(locNeigh);
                end
            end
            Hgo = mapObj.attributesNames(sum(qv(:,:,dim3) < qthresh)>0);
            H = GOind * mapObj.features;
            H = bsxfun(@rdivide, H, sum(GOind,2));
            hfig = clusterMap(H(:,sum(abs(H)>=3,1)>0),...
                Hgo,mapObj.featureNames(sum(abs(H)>=3,1)>0),3);
%             hfig = clustergram(H,'RowLabels',Hgo,'ColumnLabels',...
%                 mapObj.featureNames,'DisplayRange',6);
%             hfig.Colormap = flipud(othercolor('RdBu9'));
        end
        
        function contours = getContours(mapObj)
            ind = mapObj.contours;
            contours = ind;
            Xbins = getXbins(mapObj);
            Ybins = getYbins(mapObj);
            for kk=1:2
                for ii=1:size(ind,1)
                    for jj=1:size(ind{ii,kk},1)
                        contours{ii,kk}{jj} = [Xbins(ind{ii,kk}{jj}(:,1)),...
                            Ybins(ind{ii,kk}{jj}(:,2))];
                    end
                end
            end
        end

        
    end % End Public methods
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %       PRIVATE CLASS METHODS
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods %(Access = private)
        
        function mapObj = featureEnrichement(mapObj)
            nb = size(mapObj.features,2);
            mapObj.pvalFeatures = zeros(length(mapObj.map), nb, 2);
            % Enrichments for continuous atributes
            % Sum scores within the neighborhoods
            neigh = getNeighbors(mapObj);
            Scores = neigh * mapObj.features;
            % compute empirical backgound scores
            ScRnd = zeros(size(Scores,1), size(Scores,2), 1000);
            w = waitbar(0,{'Background features scores';'Bootstrap features scores'});
            for ii = 1:1000
                ScRnd(:,:,ii) = neigh * mapObj.features(randperm(length(mapObj.map))',:);
                waitbar(ii/1000);
            end
            waitbar(0, w, 'Calculating features background density(ies)');
            for jj=1:size(mapObj.features,2)
                ScRnd_mean = nanmean(ScRnd(:,jj,:),3);
                ScRnd_std = nanstd(ScRnd(:,jj,:),[],3);
                Z = (Scores(:,jj) - ScRnd_mean)./ScRnd_std;
                mapObj.pvalFeatures(:,jj,1) = 1-normcdf(Z); % 'high'
                mapObj.pvalFeatures(:,jj,2) = normcdf(Z);   % 'low'
                waitbar(ii/size(mapObj.features,2));
            end
            close(w);

        end
            
        function mapObj = attributeEnrichement(mapObj)

            if isempty(mapObj.attributesStyle)
                mapObj.attributesStyle = 'binary';
                U = unique(mapObj.attributes(:));
                if length(U)==3
                    mapObj.attributesStyle = 'ternary';
                elseif length(U)>3
                    mapObj.attributesStyle = 'continuous';
                end
            end

            switch mapObj.attributesStyle
                case 'binary'
                    nb = size(mapObj.attributes,2);
                    mapObj.pvalAttributes = zeros(length(mapObj.map),nb,2);
                    % enrichments for binary attributes
                    No = zeros(length(mapObj.map), size(mapObj.attributes,2)) + length(mapObj.map);
                    Ng = repmat(sum(mapObj.attributes), length(mapObj.map),1);
                    Ni = repmat(sum(mapObj.neighbors,2), 1, size(mapObj.attributes,2));
                    Nig = mapObj.neighbors * mapObj.attributes;
            %         mapObj.folds = (Nig./Ni) ./ (Ng./No);

                    w = waitbar(0,'Calculating attribute(s) enrichment(s)');
                    for ii=1:size(mapObj.attributes,2)
                        mapObj.pvalAttributes(:,ii,1) = 1 - hygecdf(Nig(:,ii)-1, No(:,ii), Ng(:,ii), Ni(:,ii));
                        mapObj.pvalAttributes(:,ii,2) = hygecdf(Nig(:,ii)-1, No(:,ii), Ng(:,ii), Ni(:,ii));
                        waitbar(ii/size(mapObj.attributes,2));
                    end
                    close(w);
                    
                case 'ternary' %{-1; 0; 1}
                    nb = size(mapObj.attributes,2);
                    mapObj.pvalAttributes = zeros(length(mapObj.map),nb,2);
                    % enrichments for top attributes
                    topAttributes = double(mapObj.attributes>0);
                    No = zeros(length(mapObj.map), size(topAttributes,2)) + length(mapObj.map);
                    Ng = repmat(sum(topAttributes), length(mapObj.map),1);
                    Ni = repmat(sum(mapObj.neighbors,2), 1, size(topAttributes,2));
                    Nig = mapObj.neighbors * topAttributes;
                    w = waitbar(0,'Calculating attribute(s) enrichment(s)');
                    for ii=1:size(topAttributes,2)
                        mapObj.pvalAttributes(:,ii,1) = 1 - hygecdf(Nig(:,ii)-1, No(:,ii), Ng(:,ii), Ni(:,ii));
                        waitbar(ii/(2*size(topAttributes,2)));
                    end
                    
                    % enrichments for bottom attributes
                    botAttributes = double(mapObj.attributes<0);
                    No = zeros(length(mapObj.map), size(botAttributes,2)) + length(mapObj.map);
                    Ng = repmat(sum(botAttributes), length(mapObj.map),1);
                    Ni = repmat(sum(mapObj.neighbors,2), 1, size(botAttributes,2));
                    Nig = mapObj.neighbors * botAttributes;
                    for ii=1:size(botAttributes,2)
                        mapObj.pvalAttributes(:,ii,2) = 1 - hygecdf(Nig(:,ii)-1, No(:,ii), Ng(:,ii), Ni(:,ii));
                        waitbar((ii+size(botAttributes,2))/(2*size(botAttributes,2)));
                    end
                    close(w);
                    
                case 'continuous'
                    nb = size(mapObj.attributes,2);
                    mapObj.pvalAttributes = zeros(length(mapObj.map), nb, 2);
                    % Enrichments for continuous atributes
                    % Sum scores within the neighborhoods
                    Scores = mapObj.neighbors * mapObj.attributes;
                    % compute empirical backgound scores
                    ScRnd = zeros(size(Scores,1), size(Scores,2), 1000);
                    w = waitbar(0,{'Background attribute(s) score(s)';'Bootstrap attribute(s) scores'});
                    for ii = 1:1000
                        ScRnd(:,:,ii) = mapObj.neighbors * mapObj.attributes(randperm(length(mapObj.map))',:);
                        waitbar(ii/1000);
                    end
                    waitbar(0, w, 'Calculating attribute(s) background density(ies)');
                    for jj=1:size(mapObj.attributes,2)
                        ScRnd_mean = nanmean(ScRnd(:,jj,:),3);
                        ScRnd_std = nanstd(ScRnd(:,jj,:),[],3);
                        Z = (Scores(:,jj) - ScRnd_mean)./ScRnd_std;
                        mapObj.pvalAttributes(:,jj,1) = 1-normcdf(Z); % 'high'
                        mapObj.pvalAttributes(:,jj,2) = normcdf(Z);   % 'low'
                        waitbar(ii/size(mapObj.attributes,2));
                    end
                    close(w);
                otherwise
                    error('Invalid definition for property attributesStyle: {continuous | binary}');
            end



        end
        
        function mapObj = correctFDRfeatures(mapObj, varargin)

            fdr = 0.05;
            if nargin ==2
                fdr = varargin{1};
            end

            pval = getPval(mapObj,'features');

            mapObj.qvalFeatures = zeros(size(pval,1),size(pval,2), size(pval,3));

            w = waitbar(0,sprintf('Correcting features enrichment p-values: FDR %2.2f%',fdr));
            for jj=1:2
            for ii = 1:size(mapObj.features,2)
                [~, mapObj.qvalFeatures(:,ii,jj)] = fdrbh(pval(:,ii,jj), fdr, 'dep');
                waitbar((jj-1)/2 + ii/(2*size(mapObj.features,2)));
            end
            end
            close(w);
            mapObj.enrNeighFeatures = mapObj.qvalFeatures <= fdr;

        end

        function mapObj = correctFDRattributes(mapObj, varargin)

            fdr = 0.05;
            if nargin ==2
                fdr = varargin{1};
            end

            pval = getPval(mapObj,'attributes');


            mapObj.qvalAttributes = zeros(size(pval,1), size(pval,2), size(pval,3));

            w = waitbar(0,sprintf('Correcting attributes enrichment p-values: FDR %2.2f%',fdr));
            for jj=1:2
            for ii = 1:size(mapObj.attributes,2)
                [~, mapObj.qvalAttributes(:,ii,jj)] = fdrbh(pval(:,ii,jj), fdr, 'dep');
                waitbar((jj-1)/2 +ii/(2*size(mapObj.attributes,2)));
            end
            end
            close(w);
            mapObj.enrNeighAttributes = mapObj.qvalAttributes <= fdr;
            fprintf('%d attributes with enriched neighborhood\n',...
                sum(sum(mapObj.enrNeighAttributes(:,:,1))>0));
            if strcmpi('continuous',mapObj.attributesStyle)
                fprintf('%d attributes with empoverished neighborhood\n',...
                    sum(sum(mapObj.enrNeighAttributes(:,:,2))>0));
            end
        end
        
        function mapObj = getFeatureContours(mapObj, varargin) 

            sigThresh = 0.05;
            if nargin == 2
                sigThresh = varargin{1};
            end
            Xbins = getXbins(mapObj);
            Ybins = getYbins(mapObj);
            qval = getQval(mapObj, 'features');
            selectNeigh = qval <= sigThresh;

            se = strel('disk', round(mapObj.radius));
            mapObj.contours = cell(size(mapObj.features,2), 2);

            for jj=1:2
                nn = find(sum(selectNeigh(:,:,jj))>0);
                for ii = nn
                    M = hist3(mapObj.map(selectNeigh(:,ii,jj),:), 'edges',...
                        {Xbins, Ybins});
                    M = imdilate(M, se);
                    set(gca,'ydir','normal');
                    msk = double(M>0);
                    mapObj.contours{ii,jj} = bwboundaries(msk, 'nohole');
                end
            end 
            close(gcf);
        end
        
        function Xbins = getXbins(mapObj)
            [~,~,X,~]=kde2dv2(mapObj.map,2^8);
            Xbins = unique(X(:));
        end
        
        function Ybins = getYbins(mapObj)
            [~,~,~,Y]=kde2dv2(mapObj.map,2^8);
            Ybins = unique(Y(:));
        end
        
    end % End Private methods
end % Class object