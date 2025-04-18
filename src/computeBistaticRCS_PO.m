function sigma_bi = computeBistaticRCS_PO(k, r, h, theta, phi1, phi2)
    % Handle the sine term fraction:
    denom = k * h * (sin(phi1) + sin(phi2));
    if abs(denom) < 1e-12
        fractionTerm = 1;
    else
        fractionTerm = sin(denom) / (denom);
    end

    % Calculate the bracket term from the paper.
    bracketTerm = ( (cos(phi2)^2 * cos(theta/2)^2) / (max(cos(phi1)^2, 1e-12)) ) * fractionTerm;
    
    % Square the entire bracket term.
    sigma_bi = k * r * h^2 * (bracketTerm)^2;
end
