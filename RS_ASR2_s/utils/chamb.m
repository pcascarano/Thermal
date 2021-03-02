function [Inew,xi]=chamb(xi,I,lambda,niter1,tau)

for j=1:niter1
        % Chambolle's iteration
        gdv = grad(div(xi) - I/lambda ); %modified to work with namespace +chambolle
        %gdv=chambolle.grad(filters.div(xi))-grad(I/lambda);
        d = sqrt(sum(gdv.^2,3));
        d = repmat( d, [1 1 2] );
        xi = ( xi + tau*gdv ) ./ ( 1+tau*d );
        % reconstruction
        
end
  Inew = I - lambda*div( xi );