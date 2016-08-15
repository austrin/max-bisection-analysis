function v = random_configuration;

A =  [-1    -1    -1
      -1     1    1
       1    -1    1
       1     1    -1];
B = ones(4, 1);

v = [cos(pi*rand(3,1))];
while (min(A*v < B) == 0)
   v = [cos(pi*rand(3,1))];
end
