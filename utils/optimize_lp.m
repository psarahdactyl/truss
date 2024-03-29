function [x,ar,ax,be,fval] = optimize_lp(obj,A,b,Aeq,beq,solver,varargin)
  % OPTIMIZE_LP 
  % 
  %
  % Inputs:
  %   (m = #E)
  %   obj 3*m by 1 vector containing the objective function value for each
  %     rod (options are: rod lengths, rod lengths * f(rod visibilities)
  %   A   4*m by 3*m inequality matrix (stress bounds)
  %   b   4*m by 1 vector of 0's
  %   Aeq #V*dim by 3*m equality matrix (force balance)
  %   beq #V*dim by 1 vectorized force matrix
  %   fval  objective of output
  

  % Outputs:
  %   (m = #E)
  %   x  3*m by 1 vector of solution of decision variables
  %   ar m by 1 vector of rod areas (1:m of x)
  %   ax m by 1 vector of rod axial forces (m+1:2*m of x)
  %   be m by 1 vector of rod bending forces (2*m+1:3*m of x)

   bending = 1;
  
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Bending'}, ...
    {'bending'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end
  
  m = size(obj,1);
  nf = 1;
  
 % set up objective
 if(bending==1)
  f = [obj;zeros(3*m*nf,1)];
 elseif(bending==0)
  f = [obj;zeros(m*nf,1)];
 end
 
 % set upper and lower bounds
 lb = [zeros(m,1);-inf(2*m,1)];
 ub = [inf*ones(m,1);inf(2*m,1)];
 
 switch solver
     case 'yalmip'
       if(bending==1)
         xo = sdpvar(4*m,1);
         Constraints = [A * xo <= b, Aeq * xo == beq];
         Objective = f' * xo;
         options = sdpsettings('verbose',1);
         sol = optimize(Constraints,Objective,options);
         save('everything.mat','A','b','Aeq','beq','f');
         
         % Analyze error flags
         if sol.problem == 0
           % Extract and display value
           x = value(xo);
           ar = x(1:m);
           ax = reshape(x(m+(1:m*nf)),m,nf);
           be = reshape(x(2*m+(1:2*m*nf)),2*m,nf);
         else
           display('Hmm, something went wrong!');
           sol.info
           yalmiperror(sol.problem)
         end
       elseif(bending==0)
         xo = sdpvar(2*m,1);
         Constraints = [A * xo <= b, Aeq * xo == beq];
         Objective = f' * xo;
         options = sdpsettings('verbose',1);

         sol = optimize(Constraints,Objective,options);

         % Analyze error flags
         if sol.problem == 0
           % Extract and display value
           x = value(xo);
           ar = x(1:m);
           ax = reshape(x(m+(1:m*nf)),m,nf);
           be =[];
         else
           display('Hmm, something went wrong!');
           sol.info
           yalmiperror(sol.problem)
         end
       end

     case 'cvx'
        cvx_begin
%         cvx_solver mosek
          variable x(4*m);
          minimize( f' * x );
          subject to
            A * x <= b;
            Aeq * x == beq;
        cvx_end
        ar = x(1:m);
        ax = reshape(x(m+(1:m*nf)),m,nf);
        be = reshape(x(2*m+(1:2*m*nf)),2*m,nf);
        
     case 'linprog'
       tic;
        [x,fval,flags,output] = linprog(f,A,b,...
            Aeq,beq);%, ...
%             lb,ub);
      fprintf('linprog: %g secs\n',toc);
        if(~isempty(x))
          ar = x(1:m);
          ax = reshape(x(m+(1:m*nf)),m,nf);
          if(bending==1)
            be = reshape(x(2*m+(1:2*m*nf)),m,2,nf);
          else
            be=[]
          end
        else
          x = [];
          ar = [];
          ax = [];
          be = [];
        end
 end
fprintf('#rods:  %d\n',sum(max(ar,0)>1e-7 & ax>0));
fprintf('#wires: %d\n',sum(max(ar,0)>1e-7 & ax<0));

end
