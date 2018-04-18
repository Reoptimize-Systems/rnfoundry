
theta = rand () * tau();

abc_mag = 1;

abc = [abc_mag*sin(0); abc_mag*sin(tau()/3); abc_mag*sin(2*tau()/3)];

dq0_1 = abc2dq0 (abc, theta, true);

alphabeta_1 = abc2alphabeta (abc);

dq0_2 = alphabeta2dq0 (alphabeta_1, theta);

alphabeta_2 = [dq02alphabeta( dq0_1, theta ); 0];

alphabeta_3 = [dq02alphabeta( dq0_2, theta ); 0];

abc_1 = dq02abc (dq0_1, theta, true);
abc_2 = dq02abc (dq0_2, theta, true);

abc_3 = alphabeta2abc (alphabeta_1);
abc_4 = alphabeta2abc (alphabeta_2);
abc_5 = alphabeta2abc (alphabeta_3);

rowheadings = { 'abc', ...
                'abc_1', ...
                'abc_2', ...
                'abc_3', ...
                'abc_4', ...
                'abc_5', ...
                'dq0_1', ...
                'dq0_2', ...
                'alphabeta_1', ...
                'alphabeta_2', ...
                'alphabeta_3' };

data = [ abc'; 
         abc_1'; 
         abc_2';
         abc_3';
         abc_4';
         abc_5';
         dq0_1';
         dq0_2';
         alphabeta_1';
         alphabeta_2';
         alphabeta_3'; ];

displaytable ( data, {}, 20, [], rowheadings);