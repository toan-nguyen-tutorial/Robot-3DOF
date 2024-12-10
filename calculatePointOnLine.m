%% Draw Point A in Line P2P3 for tool
function A = calculatePointOnLine(P2, P3, distance)
   
    direction = P3 - P2;
    direction = direction / norm(direction); 
    A = P3 - direction * distance;
  