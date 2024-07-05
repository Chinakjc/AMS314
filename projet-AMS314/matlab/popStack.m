function [ielem, stack] = popStack(stack)
% POPSTACK : pop the last pushed element on the stack, 
% return 0 if the stack is empty
ielem = 0;
if ( length(stack) ~= 0 )
  ielem = stack(end);
  stack = stack(1:end-1); 
end