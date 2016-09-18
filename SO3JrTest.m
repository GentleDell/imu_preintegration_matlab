a = [1 1 1]';
delta = [0.5, 0.5, 0.]';
for i = 1:10
    err = norm(SO3.exp(a+delta) - SO3.exp(a)*SO3.exp(SO3.Dexp(a)*delta));
%     err = norm(SO3.log(SO3.exp(a+delta) - SO3.exp(a)*SO3.exp(delta)));SO3.log(
    fprintf('delta %f, linearization error adding Jr %f\n',delta(1),err);
    delta = delta / 5;
end
fprintf('\n');
delta = [0.5, 0.5, 0.]';
for i = 1:10
%     err = norm(SO3.log(SO3.exp(a+delta) - SO3.exp(a)*SO3.exp(SO3.Dexp(a)*delta)));SO3.log(
    err = norm(SO3.exp(a+delta) - SO3.exp(a)*SO3.exp(delta));
    fprintf('delta %f, linearization error without Jr %f\n',delta(1),err);
    delta = delta / 5;
end

