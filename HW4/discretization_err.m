function err = discretization_err(fh, frh, r, p)

err = (frh - fh) / (r^p - 1);

end
