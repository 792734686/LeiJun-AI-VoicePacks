function [ P,W ] = gauspw( n )
% function to set up gauss integration point and its weight coefficient
% n --- number of gauss integration point
% P(n) --- coordinate of gauss integration point
% W(n) --- weight coefficient
switch n
    case 2
        P(1) = -0.577350269189626;
        P(2) = 0.577350269189626;
        W(1) = 1.0;
        W(2) = 1.0;
    case 3
        P(1) = -0.774596669241483;
        P(2) = 0;
        P(3) = 0.774596669241483;
        %----------------------------
        W(1) = 0.555555555555556;
        W(2) = 0.888888888888889;
        W(3) = W(1);
    case 4
        P(1) = -0.861136311594053;
        P(2) = -0.339981043584856;
        P(3) = -P(2);
        P(4) = -P(1);
        %----------------------------
        W(1) = 0.347854845137454;
        W(2) = 0.652145154862546;
        W(3) = W(2);
        W(4) = W(1);
    case 5
        P(1) = -0.906179845938664;
        P(2) = -0.538469310105683;
        P(3) = 0.0;
        P(4) = -P(2);
        P(5) = -P(1);
        %----------------------------
        W(1) = 0.236926885056189;
        W(2) = 0.478628670499366;
        W(3) = 0.568888888888889;
        W(4) = W(2);
        W(5) = W(1);
    case 6
        P(1) = -0.932469514203151;
        P(2) = -0.661209386466265;
        P(3) = -0.238619186083197;
        P(4) = -P(3);
        P(5) = -P(2);
        P(6) = -P(1);
        %----------------------------
        W(1) = 0.171324492379170;
        W(2) = 0.360761573048139;
        W(3) = 0.467913934572691;
        W(4) = W(3);
        W(5) = W(2);
        W(6) = W(1);
    case 7
        P(1) = -0.949107912342759;
        P(2) = -0.741531185599394;
        P(3) = -0.450845151377397;
        P(4) = 0.0;
        P(5) = -P(3);
        P(6) = -P(2);
        P(7) = -P(1);
        %----------------------------
        W(1) = 0.129484966168870;
        W(2) = 0.279705391489277;
        W(3) = 0.381830050505119;
        W(4) = 0.417959183673469;
        W(5) = W(3);
        W(6) = W(2);
        W(7) = W(1);
    case 8
        P(1) = -0.9602898564975363;
        P(2) = -0.7966664774136268;
        P(3) = -0.525532409916329;
        P(4) = -0.1834346424956498;
        P(5) =  0.1834346424956498;
        P(6) =  0.525532409916329;
        P(7) =  0.7966664774136268;
        P(8) =  0.9602898564975363;
        %-------------------------------
        W(1) =  0.1012285362903768;
        W(2) =  0.2223810344533745;
        W(3) =  0.3137066458778874;
        W(4) =  0.362683783378362;
        W(5) =  0.362683783378362;
        W(6) =  0.3137066458778874;
        W(7) =  0.2223810344533745;
        W(8) =  0.1012285362903768;
    case 16
        P(1) = -0.9894009349916499;
        P(2) = -0.9445750230732326;
        P(3) = -0.8656312023878318;
        P(4) = -0.755404408355003;
        P(5) = -0.6178762444026438;
        P(6) = -0.4580167776572274;
        P(7) = -0.2816035507792589;
        P(8) = -0.09501250983763744;
        P(9) =  0.09501250983763744;
        P(10) = 0.2816035507792589;
        P(11) = 0.4580167776572274;
        P(12) = 0.6178762444026438;
        P(13) = 0.755404408355003;
        P(14) = 0.8656312023878318;
        P(15) = 0.9445750230732326;
        P(16) = 0.9894009349916499;
        %-------------------------------------------
        W(1) = 0.02715245941175406;
        W(2) = 0.06225352393864778;
        W(3) = 0.0951585116824929;
        W(4) = 0.1246289712555339;
        W(5) = 0.1495959888165768;
        W(6) = 0.1691565193950026;
        W(7) = 0.1826034150449236;
        W(8) = 0.1894506104550685;
        W(9) = 0.1894506104550685;
        W(10) = 0.1826034150449236;
        W(11) = 0.1691565193950026;
        W(12) = 0.1495959888165768;
        W(13) = 0.1246289712555339;
        W(14) = 0.0951585116824929;
        W(15) = 0.06225352393864778;
        W(16) = 0.02715245941175406;
end
end
