function [GPNatCoords, GPweights] = gaussPoints(elementType,numGP)
% Coordinates and weights of Gauss points for the specified element type
% and number of points.
%
% INPUT:
% elementType - a string specifying the type of element.
% numGP - single scalar specifying the requested number of Gauss points.
%
% OUTPUT:
% GPNatCoords - #GP-by-#ElDim  matrix,  row GPNatCoords(k,:) contains the
%   coordinates of k-th Gauss point.
% GPweights - #GP-by-1 matrix, weights(k,1) is the k-th Gauss weight.
%
% Examples:
% [GPNatCoords, GPweights]=gaussPoints('2dQ1',4)
%
% COMMENTS:
% Supported element types and corresponding numbers of points:
% '1dQ1', '1dQ1': 1, 2, 3, 4 
% '1dP1', '1dP1': 1, 2, 3, 4 
% '2dQ1', '2dQ2', '2dQ2r': 1, 4, 9, 16
% '2dP1', '2dP2': 3, 4, 6, 7
% '3dQ1', '3dQ2', '3dQ2r': 1, 8, 27
% '3dP1', '3dP2': 1

if nargin<2
    numGP=1;
end

if elementType(1)=='1' %1D case
    
    switch elementType
        case {'1dQ1','1dQ2'}
            if numGP==1
                GPNatCoords=0;
                GPweights=2;
            elseif numGP==2
                GPNatCoords=[-0.577350269189626; 0.577350269189626];
                GPweights=[1 1];
            elseif numGP==3
                GPNatCoords=[-0.774596669241483; 0.774596669241483; 0];
                GPweights=[0.555555555555556; 0.555555555555556; 0.888888888888889];
            elseif numGP==4
                GPNatCoords=[0.339981043584856; -0.339981043584856; 0.861136311594053; 0.861136311594053];
                GPweights=[0.652145154862546; 0.652145154862546; 0.347854845137454; 0.347854845137454];
            else
                error('Number of points specified is not supported');
            end
        case {'1dP1','1dP2'}
            if numGP==1
                GPNatCoords=0;
                GPweights=2/2;
            elseif numGP==2
                GPNatCoords=[-0.577350269189626; 0.577350269189626];
                GPweights=[1 1]/2;
            elseif numGP==3
                GPNatCoords=[-0.774596669241483; 0.774596669241483; 0];
                GPweights=[0.555555555555556; 0.555555555555556; 0.888888888888889]/2;
            elseif numGP==4
                GPNatCoords=[0.339981043584856; -0.339981043584856; 0.861136311594053; 0.861136311594053];
                GPweights=[0.652145154862546; 0.652145154862546; 0.347854845137454; 0.347854845137454]/2;
            else
                error('Number of points specified is not supported');
            end
    end
    
elseif elementType(1)=='2' %2D case

    switch elementType
        case {'2dQ1','2dQ2','2dQ2r'}
            if numGP==1
                GPNatCoords=[0.000000000000000e+00, 0.000000000000000e+00];
                GPweights=[4.000000000000000e+00];
            elseif numGP==4
                GPNatCoords=[-5.773502691896257e-01, -5.773502691896257e-01;
                    5.773502691896257e-01, -5.773502691896257e-01;
                    -5.773502691896257e-01, 5.773502691896257e-01;
                    5.773502691896257e-01, 5.773502691896257e-01];
                GPweights=[1.000000000000000e+00;
                    1.000000000000000e+00;
                    1.000000000000000e+00;
                    1.000000000000000e+00];
            elseif numGP==9
                GPNatCoords=[0.000000000000000e+00, 0.000000000000000e+00;
                    -7.745966692414834e-01, 0.000000000000000e+00;
                    7.745966692414834e-01, 0.000000000000000e+00;
                    0.000000000000000e+00, -7.745966692414834e-01;
                    -7.745966692414834e-01, -7.745966692414834e-01;
                    7.745966692414834e-01, -7.745966692414834e-01;
                    0.000000000000000e+00, 7.745966692414834e-01;
                    -7.745966692414834e-01, 7.745966692414834e-01;
                    7.745966692414834e-01, 7.745966692414834e-01];
                GPweights=[7.901234567901234e-01;
                    4.938271604938271e-01;
                    4.938271604938271e-01;
                    4.938271604938271e-01;
                    3.086419753086420e-01;
                    3.086419753086420e-01;
                    4.938271604938271e-01;
                    3.086419753086420e-01;
                    3.086419753086420e-01];
            elseif numGP==16
                GPNatCoords=[-3.399810435848563e-01, -3.399810435848563e-01;
                    3.399810435848563e-01, -3.399810435848563e-01;
                    -8.611363115940526e-01, -3.399810435848563e-01;
                    8.611363115940526e-01, -3.399810435848563e-01;
                    -3.399810435848563e-01, 3.399810435848563e-01;
                    3.399810435848563e-01, 3.399810435848563e-01;
                    -8.611363115940526e-01, 3.399810435848563e-01;
                    8.611363115940526e-01, 3.399810435848563e-01;
                    -3.399810435848563e-01, -8.611363115940526e-01;
                    3.399810435848563e-01, -8.611363115940526e-01;
                    -8.611363115940526e-01, -8.611363115940526e-01;
                    8.611363115940526e-01, -8.611363115940526e-01;
                    -3.399810435848563e-01, 8.611363115940526e-01;
                    3.399810435848563e-01, 8.611363115940526e-01;
                    -8.611363115940526e-01, 8.611363115940526e-01;
                    8.611363115940526e-01, 8.611363115940526e-01];
                GPweights=[4.252933030106944e-01;
                    4.252933030106944e-01;
                    2.268518518518519e-01;
                    -8.110558238879494e-02;
                    4.252933030106944e-01;
                    4.252933030106944e-01;
                    2.268518518518519e-01;
                    -8.110558238879494e-02;
                    2.268518518518519e-01;
                    2.268518518518519e-01;
                    1.210029932856020e-01;
                    -4.326179469597342e-02;
                    -8.110558238879494e-02;
                    -8.110558238879494e-02;
                    -4.326179469597342e-02;
                    1.546724448294497e-02];
            else
                error('Number of points specified is not supported');
            end
        case {'2dP1','2dP2'}
            if numGP==1
                GPNatCoords=[0.3333333333333 0.3333333333333];
                GPweights=0.5;
            elseif numGP==3
                GPNatCoords=[0.1666666666667, 0.1666666666667;
                    0.6666666666667, 0.1666666666667;
                    0.1666666666667, 0.6666666666667];
                GPweights=[0.1666666666667;0.1666666666667;0.1666666666667];
            elseif numGP==4
                GPNatCoords=[.6, .2;
                    .2, .6;
                    .2, .2;
                    .3333333333333, .3333333333333];
                GPweights=[0.260416666666667
                    0.260416666666667
                    0.260416666666667
                    -0.281250000000000];
            elseif numGP==6
                GPNatCoords=[...
                    0.816847572980459, 0.091576213509771;
                    0.091576213509771, 0.816847572980459;
                    0.091576213509771, 0.091576213509771;
                    0.108103018168070, 0.445948490915965;
                    0.445948490915965, 0.108103018168070;
                    0.445948490915965, 0.445948490915965];
                GPweights=[...
                    0.109951743655322;
                    0.109951743655322;
                    0.109951743655322;
                    0.223381589678011;
                    0.223381589678011;
                    0.223381589678011];

            elseif numGP==7
                GPNatCoords=[0.1012865073235, 0.1012865073235;
                    0.7974269853531, 0.1012865073235;
                    0.1012865073235, 0.7974269853531;
                    0.4701420641051, 0.0597158717898;
                    0.4701420641051, 0.4701420641051;
                    0.0597158717898, 0.4701420641051;
                    0.3333333333333, 0.3333333333333];
                GPweights=[0.062969590272400;
                    0.062969590272400;
                    0.062969590272400;
                    0.066197076394250;
                    0.066197076394250;
                    0.066197076394250;
                    0.112500000000000];
            else
                error('Number of points specified is not supported for the chosen type of element');
            end
        otherwise
            error('Element type is not supported')
    end
    
else %3D case
    
    switch elementType
        case {'3dQ1','3dQ2','3dQ2r'}
            if numGP==1
                GPNatCoords=[0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00];
                GPweights=[8.000000000000000e+00];
            elseif numGP==8
                GPNatCoords=[-5.773502691896257e-01, -5.773502691896257e-01, -5.773502691896257e-01;
                    5.773502691896257e-01, -5.773502691896257e-01, -5.773502691896257e-01;
                    -5.773502691896257e-01, 5.773502691896257e-01, -5.773502691896257e-01;
                    5.773502691896257e-01,  5.773502691896257e-01, -5.773502691896257e-01;
                    -5.773502691896257e-01, -5.773502691896257e-01, 5.773502691896257e-01;
                    5.773502691896257e-01, -5.773502691896257e-01, 5.773502691896257e-01;
                    -5.773502691896257e-01, 5.773502691896257e-01, 5.773502691896257e-01;
                    5.773502691896257e-01, 5.773502691896257e-01, 5.773502691896257e-01];
                GPweights=[1.000000000000000e+00;
                    1.000000000000000e+00;
                    1.000000000000000e+00;
                    1.000000000000000e+00;
                    1.000000000000000e+00;
                    1.000000000000000e+00;
                    1.000000000000000e+00;
                    1.000000000000000e+00];

            elseif numGP==27
                GPNatCoords=[0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00;
                    -7.745966692414834e-01, 0.000000000000000e+00, 0.000000000000000e+00;
                    7.745966692414834e-01, 0.000000000000000e+00, 0.000000000000000e+00;
                    0.000000000000000e+00, -7.745966692414834e-01, 0.000000000000000e+00;
                    -7.745966692414834e-01, -7.745966692414834e-01, 0.000000000000000e+00;
                    7.745966692414834e-01, -7.745966692414834e-01, 0.000000000000000e+00;
                    0.000000000000000e+00, 7.745966692414834e-01, 0.000000000000000e+00;
                    -7.745966692414834e-01, 7.745966692414834e-01, 0.000000000000000e+00;
                    7.745966692414834e-01, 7.745966692414834e-01, 0.000000000000000e+00;
                    0.000000000000000e+00, 0.000000000000000e+00, -7.745966692414834e-01;
                    -7.745966692414834e-01, 0.000000000000000e+00, -7.745966692414834e-01;
                    7.745966692414834e-01, 0.000000000000000e+00, -7.745966692414834e-01;
                    0.000000000000000e+00, -7.745966692414834e-01, -7.745966692414834e-01;
                    -7.745966692414834e-01, -7.745966692414834e-01, -7.745966692414834e-01;
                    7.745966692414834e-01, -7.745966692414834e-01, -7.745966692414834e-01;
                    0.000000000000000e+00, 7.745966692414834e-01, -7.745966692414834e-01;
                    -7.745966692414834e-01, 7.745966692414834e-01, -7.745966692414834e-01;
                    7.745966692414834e-01, 7.745966692414834e-01, -7.745966692414834e-01;
                    0.000000000000000e+00, 0.000000000000000e+00, 7.745966692414834e-01;
                    -7.745966692414834e-01, 0.000000000000000e+00, 7.745966692414834e-01;
                    7.745966692414834e-01, 0.000000000000000e+00, 7.745966692414834e-01;
                    0.000000000000000e+00, -7.745966692414834e-01, 7.745966692414834e-01;
                    -7.745966692414834e-01, -7.745966692414834e-01, 7.745966692414834e-01;
                    7.745966692414834e-01, -7.745966692414834e-01, 7.745966692414834e-01;
                    0.000000000000000e+00, 7.745966692414834e-01, 7.745966692414834e-01;
                    -7.745966692414834e-01, 7.745966692414834e-01, 7.745966692414834e-01;
                    7.745966692414834e-01, 7.745966692414834e-01, 7.745966692414834e-01];
                GPweights=[7.023319615912208e-01;4.389574759945130e-01;
                    4.389574759945130e-01;4.389574759945130e-01;
                    2.743484224965707e-01;2.743484224965707e-01;
                    4.389574759945130e-01;2.743484224965707e-01;
                    2.743484224965707e-01;4.389574759945130e-01;
                    2.743484224965707e-01;2.743484224965707e-01;
                    2.743484224965707e-01;1.714677640603567e-01;
                    1.714677640603567e-01;2.743484224965707e-01;
                    1.714677640603567e-01;1.714677640603567e-01;
                    4.389574759945130e-01;2.743484224965707e-01;
                    2.743484224965707e-01;2.743484224965707e-01;
                    1.714677640603567e-01;1.714677640603567e-01;
                    2.743484224965707e-01;1.714677640603567e-01;
                    1.714677640603567e-01];
            end
        case {'3dP1','3dP2'}
            if numGP==1
                GPNatCoords=0.25*[1 1 1];
                GPweights=[1/6];
            else
                error('Number of points specified is not supported');
            end
        otherwise
            error('Element type is not supported')
    end
    
end
