sor_2d_hom = '''
	# line 13 "sor_poisson_2d_hom.py"
	double tmp, diff, err=1.0, sum, dx2 = double(_dx2);
	int t=0, i, j, k;

	for (t=0; t<int(nmax); t++) {
		err = 0.0;
		sum = 0.0;

	    i=0, j=0, k=0; // left front bottom corner neumann boundary conditions with 1st order approximation
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*u(i,j+1,k) + I(i,j,k)*dx2)/(sxx + syy);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		i=0, k=0; // left bottom edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*(u(i,j+1,k) + u(i,j-1,k)) + I(i,j,k)*dx2)/(sxx + 2*syy);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=0, j=ny-1, k=0; // left back bottom corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*u(i,j-1,k) + I(i,j,k)*dx2)/(sxx + syy);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		j=0, k=0; // front bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*u(i,j+1,k) + I(i,j,k)*dx2)/(2*sxx + syy);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		k=0; // bottom face
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*(u(i,j+1,k) + u(i,j-1,k)) + I(i,j,k)*dx2)/(2*sxx + 2*syy);
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		j=ny-1, k=0; // back bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*u(i,j-1,k) + I(i,j,k)*dx2)/(2*sxx + syy);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, j=0, k=0; // bottom front right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(syy*u(i,j+1,k) + sxx*u(i-1,j,k) + I(i,j,k)*dx2)/(syy + sxx);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		i=nx-1, k=0; // bottom right edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(syy*(u(i,j+1,k) + u(i,j-1,k)) + sxx*u(i-1,j,k) + I(i,j,k)*dx2)/(2*syy + sxx);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}
		
		i=nx-1, j=ny-1, k=0; // bottom back right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i-1,j,k) + syy*u(i,j-1,k) + I(i,j,k)*dx2)/(sxx + syy);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		if (fabs(sum) != 0.) {
	        err = sqrt(err/sum);
	        if (err<double(tol)) break;
	    }  
	}

	if (err>double(tol)) return_val=double(err);
	else return_val = double(t);
	'''

sor_2d_hom_xy = '''
	# line 13 "sor_poisson_2d_hom_xy.py"
	double tmp, diff, err=1.0, sum, dx2 = double(_dx2);
	int t=0, i, j, k;

	for (t=0; t<int(nmax); t++) {
		err = 0.0;
		sum = 0.0;

	    i=0, j=0, k=0; // left front bottom corner neumann boundary conditions with 1st order approximation
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+sxx*u(i+1,j,k)
									   +syy*u(i,j+1,k)
									   +0.5*sxy*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
									   -0.5*sxy*(u(i,j+1,k)-u(i+1,j+1,k)-u(i+1,j,k))
									   +I(i,j,k)*dx2)/(
									   +sxx + sxy
									   +syy);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		i=0, k=0; // left bottom edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+sxx*u(i+1,j,k)
										   +syy*u(i,j+1,k)
										   +syy*u(i,j-1,k)
										   +0.25*sxy*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))
										   -0.5*sxy*(u(i,j+1,k)-u(i+1,j+1,k)-u(i+1,j,k))
										   -0.5*sxy*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   +sxx
										   +syy
										   +syy);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=0, j=ny-1, k=0; // left back bottom corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+sxx*u(i+1,j,k)
									   +syy*u(i,j-1,k)
									   +0.5*sxy*(u(i+1,j,k)-u(i,j-1,k)-u(i+1,j-1,k))
									   -0.5*sxy*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j,k))
									   +I(i,j,k)*dx2)/(
									   +sxx - sxy
									   +syy);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		j=0, k=0; // front bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+sxx*u(i+1,j,k)
										   +syy*u(i,j+1,k)
										   +sxx*u(i-1,j,k)
										   +0.5*sxy*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
										   -0.25*sxy*(u(i-1,j+1,k)+u(i-1,j,k)-u(i+1,j+1,k)-u(i+1,j,k))
										   +0.5*sxy*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   +I(i,j,k)*dx2)/(
										   2*sxx
										   +syy);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		k=0; // bottom face
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+sxx*(u(i+1,j,k) + u(i-1,j,k))
											   +syy*(u(i,j+1,k) + u(i,j-1,k))
											   +0.5*sxy*(u(i+1,j+1,k) + u(i-1,j-1,k) -u(i+1,j-1,k) - u(i-1,j+1,k))
											   +I(i,j,k)*dx2)/(
											   2*sxx + 2*syy);
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		j=ny-1, k=0; // back bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+sxx*u(i+1,j,k)
										   +sxx*u(i-1,j,k)
										   +syy*u(i,j-1,k)
										   +0.5*sxy*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
										   +0.5*sxy*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   -0.25*sxy*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j-1,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   2*sxx + syy);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, j=0, k=0; // bottom front right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+syy*u(i,j+1,k)
									   +sxx*u(i-1,j,k)
									   -0.5*sxy*(u(i-1,j+1,k)+u(i-1,j,k)-u(i,j+1,k))
									   +0.5*sxy*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
									   +I(i,j,k)*dx2)/(
									   syy - sxy + sxx);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		i=nx-1, k=0; // bottom right edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+syy*u(i,j+1,k)
										   +sxx*u(i-1,j,k)
										   +syy*u(i,j-1,k)
										   -0.5*sxy*(u(i-1,j+1,k)+u(i-1,j,k)-u(i,j+1,k))
										   +0.25*sxy*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   -0.5*sxy*(u(i,j-1,k)-u(i-1,j-1,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   2*syy + sxx);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}
		
		i=nx-1, j=ny-1, k=0; // bottom back right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+sxx*u(i-1,j,k)
									   +syy*u(i,j-1,k)
									   +0.5*sxy*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j,k))
									   -0.5*sxy*(u(i,j-1,k)-u(i-1,j-1,k)-u(i-1,j,k))
									   +I(i,j,k)*dx2)/(
									   sxx + sxy + syy);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		if (fabs(sum) != 0.) {
	        err = sqrt(err/sum);
	        if (err<double(tol)) break;
	    }  
	}

	if (err>double(tol)) return_val=double(err);
	else return_val = double(t);
	'''

sor_2d_inh_xy = '''
	# line 13 "sor_poisson_2d_inh_xy.py"
	double tmp, diff, err=1.0, sum, dx2 = double(_dx2);
	int t=0, i, j, k;

	for (t=0; t<int(nmax); t++) {
		err = 0.0;
		sum = 0.0;

	    i=0, j=0, k=0; // left front bottom corner neumann boundary conditions with 1st order approximation
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
									   +0.5*(syy(i,j+1,k)+syy(i+1,j+1,k))*u(i,j+1,k)
									   +0.25*(sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
									   -0.25*(sxy(i,j+1,k)+sxy(i+1,j+1,k))*(u(i,j+1,k)-u(i+1,j+1,k)-u(i+1,j,k))
									   +I(i,j,k)*dx2)/(
									   +0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k)) + 0.25*(sxy(i+1,j+1,k)+sxy(i+1,j,k))
									   +0.5*(syy(i,j+1,k)+syy(i+1,j+1,k)) + 0.25*(sxy(i,j+1,k)+sxy(i+1,j+1,k)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=0, k=0; // left bottom edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.5*(syy(i,j+1,k)+syy(i+1,j+1,k))*u(i,j+1,k)
										   +0.5*(syy(i,j,k)+syy(i+1,j,k))*u(i,j-1,k)
										   +0.125*(sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))
										   -0.25*(sxy(i,j+1,k)+sxy(i+1,j+1,k))*(u(i,j+1,k)-u(i+1,j+1,k)-u(i+1,j,k))
										   -0.25*(sxy(i,j,k)+sxy(i+1,j,k))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   +0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k))
										   +0.5*(syy(i,j+1,k)+syy(i+1,j+1,k)) + 0.25*(sxy(i,j+1,k)+sxy(i+1,j+1,k))
										   +0.5*(syy(i,j,k)+syy(i+1,j,k)) - 0.25*(sxy(i,j,k)+sxy(i+1,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=0, j=ny-1, k=0; // left back bottom corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
									   +0.5*(syy(i,j,k)+syy(i+1,j,k))*u(i,j-1,k)
									   +0.25*(sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i+1,j,k)-u(i,j-1,k)-u(i+1,j-1,k))
									   -0.25*(sxy(i,j,k)+sxy(i+1,j,k))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j,k))
									   +I(i,j,k)*dx2)/(
									   +0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k)) - 0.25*(sxy(i+1,j+1,k)+sxy(i+1,j,k))
									   +0.5*(syy(i,j,k)+syy(i+1,j,k)) - 0.25*(sxy(i,j,k)+sxy(i+1,j,k)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		j=0, k=0; // front bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.5*(syy(i,j+1,k)+syy(i+1,j+1,k))*u(i,j+1,k)
										   +0.5*(sxx(i,j+1,k)+sxx(i,j,k))*u(i-1,j,k)
										   +0.25*(sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
										   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i+1,j+1,k)-u(i+1,j,k))
										   +0.25*(sxy(i,j+1,k)+sxy(i,j,k))*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   +I(i,j,k)*dx2)/(
										   +0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k)) + 0.25*(sxy(i+1,j+1,k)+sxy(i+1,j,k))
										   +0.5*(syy(i,j+1,k)+syy(i+1,j+1,k))
										   +0.5*(sxx(i,j+1,k)+sxx(i,j,k)) - 0.25*(sxy(i,j+1,k)+sxy(i,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		k=0; // bottom face
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
											   +0.5*(syy(i,j+1,k)+syy(i+1,j+1,k))*u(i,j+1,k)
											   +0.5*(sxx(i,j+1,k)+sxx(i,j,k))*u(i-1,j,k)
											   +0.5*(syy(i,j,k)+syy(i+1,j,k))*u(i,j-1,k)
											   +0.125*(sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))
											   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i+1,j+1,k)-u(i+1,j,k))
											   +0.125*(sxy(i,j+1,k)+sxy(i,j,k))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j+1,k)-u(i,j+1,k))
											   -0.125*(sxy(i,j,k)+sxy(i+1,j,k))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j-1,k)-u(i-1,j,k))
											   +I(i,j,k)*dx2)/(
											   +0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k))
											   +0.5*(syy(i,j+1,k)+syy(i+1,j+1,k))
											   +0.5*(sxx(i,j+1,k)+sxx(i,j,k))
											   +0.5*(syy(i,j,k)+syy(i+1,j,k)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		j=ny-1, k=0; // back bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.5*(sxx(i,j+1,k)+sxx(i,j,k))*u(i-1,j,k)
										   +0.5*(syy(i,j,k)+syy(i+1,j,k))*u(i,j-1,k)
										   +0.25*(sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
										   +0.25*(sxy(i,j+1,k)+sxy(i,j,k))*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   -0.125*(sxy(i,j,k)+sxy(i+1,j,k))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j-1,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   +0.5*(sxx(i+1,j+1,k)+sxx(i+1,j,k)) + 0.25*(sxy(i+1,j+1,k)+sxy(i+1,j,k))
										   +0.5*(sxx(i,j+1,k)+sxx(i,j,k)) - 0.25*(sxy(i,j+1,k)+sxy(i,j,k))
										   +0.5*(syy(i,j,k)+syy(i+1,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=nx-1, j=0, k=0; // bottom front right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.5*(syy(i,j+1,k)+syy(i+1,j+1,k))*u(i,j+1,k)
									   +0.5*(sxx(i,j+1,k)+sxx(i,j,k))*u(i-1,j,k)
									   -0.25*(sxy(i,j+1,k)+sxy(i+1,j+1,k))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i,j+1,k))
									   +0.25*(sxy(i,j+1,k)+sxy(i,j,k))*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
									   +I(i,j,k)*dx2)/(
									   +0.5*(syy(i,j+1,k)+syy(i+1,j+1,k)) - 0.25*(sxy(i,j+1,k)+sxy(i+1,j+1,k))
									   +0.5*(sxx(i,j+1,k)+sxx(i,j,k)) - 0.25*(sxy(i,j+1,k)+sxy(i,j,k)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, k=0; // bottom right edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.5*(syy(i,j+1,k)+syy(i+1,j+1,k))*u(i,j+1,k)
										   +0.5*(sxx(i,j+1,k)+sxx(i,j,k))*u(i-1,j,k)
										   +0.5*(syy(i,j,k)+syy(i+1,j,k))*u(i,j-1,k)
										   -0.25*(sxy(i,j+1,k)+sxy(i+1,j+1,k))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i,j+1,k))
										   +0.125*(sxy(i,j+1,k)+sxy(i,j,k))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   -0.25*(sxy(i,j,k)+sxy(i+1,j,k))*(u(i,j-1,k)-u(i-1,j-1,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   +0.5*(syy(i,j+1,k)+syy(i+1,j+1,k)) - 0.25*(sxy(i,j+1,k)+sxy(i+1,j+1,k))
										   +0.5*(sxx(i,j+1,k)+sxx(i,j,k))
										   +0.5*(syy(i,j,k)+syy(i+1,j,k)) + 0.25*(sxy(i,j,k)+sxy(i+1,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}
		
		i=nx-1, j=ny-1, k=0; // bottom back right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.5*(sxx(i,j+1,k)+sxx(i,j,k))*u(i-1,j,k)
									   +0.5*(syy(i,j,k)+syy(i+1,j,k))*u(i,j-1,k)
									   +0.25*(sxy(i,j+1,k)+sxy(i,j,k))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j,k))
									   -0.25*(sxy(i,j,k)+sxy(i+1,j,k))*(u(i,j-1,k)-u(i-1,j-1,k)-u(i-1,j,k))
									   +I(i,j,k)*dx2)/(
									   +0.5*(sxx(i,j+1,k)+sxx(i,j,k)) + 0.25*(sxy(i,j+1,k)+sxy(i,j,k))
									   +0.5*(syy(i,j,k)+syy(i+1,j,k)) + 0.25*(sxy(i,j,k)+sxy(i+1,j,k)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		if (fabs(sum) != 0.) {
	        err = sqrt(err/sum);
	        if (err<double(tol)) break;
	    }  
	}

	if (err>double(tol)) return_val=double(err);
	else return_val = double(t);
	'''

sor_3d_hom = '''
	# line 13 "sor_poisson_3d_hom.py"
	double tmp, diff, err=1.0, sum, dx2 = double(_dx2);
	int t=0, i, j, k;

	for (t=0; t<int(nmax); t++) {
		err = 0.0;
		sum = 0.0;

	    i=0, j=0, k=0; // left front bottom corner neumann boundary conditions with 1st order approximation
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*u(i,j+1,k) + szz*u(i,j,k+1) + I(i,j,k)*dx2)/(sxx + syy + szz);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		i=0, k=0; // left bottom edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*(u(i,j+1,k) + u(i,j-1,k)) + szz*u(i,j,k+1) + I(i,j,k)*dx2)/(sxx + 2*syy + szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=0, j=ny-1, k=0; // left back bottom corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*u(i,j-1,k) + szz*(u(i,j,k+1)) + I(i,j,k)*dx2)/(
									   sxx + syy + szz);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		i=0, j=0; // front left edge
		for (k=1; k<nz-1; k++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*u(i,j+1,k) + szz*(u(i,j,k-1) + u(i,j,k+1)) + I(i,j,k)*dx2)/(
										   sxx + syy + 2*szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=0; // left face
		for (j=1; j<ny-1; j++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*(u(i,j+1,k) + u(i,j-1,k)) + szz*(u(i,j,k-1) + u(i,j,k+1)) + I(i,j,k)*dx2)/(
											   sxx + 2*syy + 2*szz);
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=0, j=ny-1; // left back edge
		for (k=1; k<nz-1; k++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*u(i,j-1,k) + szz*(u(i,j,k-1) + u(i,j,k+1)) + I(i,j,k)*dx2)/(
										   sxx + syy + 2*szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=0, j=0, k=nz-1; // left front top corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*u(i,j+1,k) + szz*u(i,j,k-1) + I(i,j,k)*dx2)/(
									   sxx + syy + szz);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		i=0, k=nz-1; // left top edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*(u(i,j+1,k) + u(i,j-1,k)) + szz*u(i,j,k-1) + I(i,j,k)*dx2)/(
										   sxx + 2*syy + szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=0, j=ny-1, k=nz-1; // left back top corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i+1,j,k) + syy*u(i,j-1,k) + szz*u(i,j,k-1) + I(i,j,k)*dx2)/(
									   sxx + syy + szz);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		j=0, k=0; // front bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*u(i,j+1,k) + szz*u(i,j,k+1) + I(i,j,k)*dx2)/(
										   2*sxx + syy + szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}
		
		int j=0; // front face
		for (i=1; i<nx-1; i++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*u(i,j+1,k) + szz*(u(i,j,k-1) + u(i,j,k+1)) + I(i,j,k)*dx2)/(
											   2*sxx + syy + 2*szz);
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		j=0,k=nz-1; // front top edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*u(i,j+1,k) + szz*u(i,j,k-1) + I(i,j,k)*dx2)/(
										   2*sxx + syy + szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		k=0; // bottom face
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*(u(i,j+1,k) + u(i,j-1,k)) + szz*u(i,j,k+1) + I(i,j,k)*dx2)/(
											   2*sxx + 2*syy +szz);
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		j=ny-1, k=0; // back bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*u(i,j-1,k) + szz*u(i,j,k+1) + I(i,j,k)*dx2)/(
										   2*sxx + syy + szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		// bulk
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				for (k=1; k<nz-1; k++) {
					tmp = u(i,j,k);
					u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*(u(i,j+1,k) + u(i,j-1,k)) + szz*(u(i,j,k-1) + u(i,j,k+1)) + I(i,j,k)*dx2)/(
												   2*sxx + 2*syy + 2*szz);
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		k=nz-1; // top face
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*(u(i,j+1,k) + u(i,j-1,k)) + szz*u(i,j,k-1) + I(i,j,k)*dx2)/(
											   2*sxx + 2*syy + szz);
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		j=ny-1; // back face
		for (i=1; i<nx-1; i++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*u(i,j-1,k) + szz*(u(i,j,k-1) + u(i,j,k+1)) + I(i,j,k)*dx2)/(
											   2*sxx + syy + 2*szz);
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		j=ny-1, k=nz-1; // back top edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*(u(i+1,j,k) + u(i-1,j,k)) + syy*u(i,j-1,k) + szz*u(i,j,k-1) + I(i,j,k)*dx2)/(
										   2*sxx + syy + szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, j=0, k=0; // bottom front right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(syy*u(i,j+1,k) + sxx*u(i-1,j,k) + szz*u(i,j,k+1) + I(i,j,k)*dx2)/(
									   syy + sxx + szz);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		i=nx-1, k=0; // bottom right edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(syy*(u(i,j+1,k) + u(i,j-1,k)) + sxx*u(i-1,j,k) + szz*u(i,j,k+1) + I(i,j,k)*dx2)/(
										   2*syy + sxx + szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}
		
		i=nx-1, j=ny-1, k=0; // bottom back right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i-1,j,k) + syy*u(i,j-1,k) + szz*u(i,j,k+1) + I(i,j,k)*dx2)/(
									   sxx + syy + szz);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		i=nx-1, j=0; // front right edge
		for (k=1; k<nz-1; k++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(syy*u(i,j+1,k) + sxx*u(i-1,j,k) + szz*(u(i,j,k-1) + u(i,j,k+1)) + I(i,j,k)*dx2)/(
										   syy + sxx + 2*szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, j=0, k=nz-1; // front right top corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(syy*u(i,j+1,k) + sxx*u(i-1,j,k) + szz*u(i,j,k-1) + I(i,j,k)*dx2)/(
									   syy + sxx + szz);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		i=nx-1; // right face
		for (j=1; j<ny-1; j++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(syy*(u(i,j+1,k) + u(i,j-1,k)) + sxx*u(i-1,j,k) + szz*(u(i,j,k-1) + u(i,j,k+1)) + I(i,j,k)*dx2)/(
											   2*syy + sxx + 2*szz);
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=nx-1, j=ny-1; // back right edge
		for (k=1; k<nz-1; k++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i-1,j,k) + syy*u(i,j-1,k) + szz*(u(i,j,k-1) + u(i,j,k+1)) + I(i,j,k)*dx2)/(
										   sxx + syy + 2*szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, k=nz-1; // right top edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(syy*(u(i,j+1,k) + u(i,j-1,k)) + sxx*u(i-1,j,k) + szz*u(i,j,k-1) + I(i,j,k)*dx2)/(
										   2*syy + sxx + szz);
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, j=ny-1, k=nz-1; // top right back corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(sxx*u(i-1,j,k) + syy*u(i,j-1,k) + szz*u(i,j,k-1) + I(i,j,k)*dx2)/(
									   sxx + syy + szz);
		diff = u(i,j,k) - tmp;
		err += diff*diff;
		sum += tmp*tmp;

		if (fabs(sum) != 0.) {
	        err = sqrt(err/sum);
	        if (err<double(tol)) break;
	    }  
	}

	if (err>double(tol)) return_val=double(err);
	else return_val = double(t);
	'''

sor_3d_inh = '''
	# line 14 "sor_poisson_3d_inh.py"
	double tmp, diff, err=1.0, sum, dx2 = double(_dx2);
	int t=0, i, j, k;

	for (t=0; t<int(nmax); t++) {
		err = 0.0;
		sum = 0.0;

		// x=0, y=0, z=0 corner neumann boundary conditions with 1st order approximation
		tmp = u(0,0,0);
		u(0,0,0) = (1-w)*u(0,0,0) + w*(+0.25*(sxx(1,1,1)+sxx(1,0,1)+sxx(1,1,0)+sxx(1,0,0))*u(1,0,0)
		                               +0.25*(syy(0,1,0)+syy(1,1,0)+syy(0,1,1)+syy(1,1,1))*u(0,1,0)
		                               +0.25*(szz(0,1,1)+szz(1,1,1)+szz(0,0,1)+szz(1,0,1))*u(0,0,1)
		                               +I(0,0,0)*dx2)/(
		                               +0.25*(sxx(1,1,1)+sxx(1,0,1)+sxx(1,1,0)+sxx(1,0,0))
		                               +0.25*(syy(0,1,0)+syy(1,1,0)+syy(0,1,1)+syy(1,1,1))
		                               +0.25*(szz(0,1,1)+szz(1,1,1)+szz(0,0,1)+szz(1,0,1)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		// x=0, z=0 edge
		for (j=1; j<ny-1; j++) {
			tmp = u(0,j,0);
			u(0,j,0) = (1-w)*u(0,j,0) + w*(+0.25*(sxx(1,j+1,1)+sxx(1,j,1)+sxx(1,j+1,0)+sxx(1,j,0))*u(1,j,0)
										   +0.25*(syy(0,j+1,0)+syy(1,j+1,0)+syy(0,j+1,1)+syy(1,j+1,1))*u(0,j+1,0)
										   +0.25*(syy(0,j,0)+syy(1,j,0)+syy(0,j,1)+syy(1,j,1))*u(0,j-1,0)
										   +0.25*(szz(0,j+1,1)+szz(1,j+1,1)+szz(0,j,1)+szz(1,j,1))*u(0,j,1)
										   +I(0,j,0)*dx2)/(
										   +0.25*(sxx(1,j+1,1)+sxx(1,j,1)+sxx(1,j+1,0)+sxx(1,j,0)) 
										   +0.25*(syy(0,j+1,0)+syy(1,j+1,0)+syy(0,j+1,1)+syy(1,j+1,1))
										   +0.25*(syy(0,j,0)+syy(1,j,0)+syy(0,j,1)+syy(1,j,1))
										   +0.25*(szz(0,j+1,1)+szz(1,j+1,1)+szz(0,j,1)+szz(1,j,1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		// x=0, y=ny, z=0 corner
		tmp = u(0,ny-1,0);
		u(0,ny-1,0) = (1-w)*u(0,ny-1,0) + w*(+0.25*(sxx(1,ny,1)+sxx(1,ny-1,1)+sxx(1,ny,0)+sxx(1,ny-1,0))*u(1,ny-1,0)
											 +0.25*(syy(0,ny-1,0)+syy(1,ny-1,0)+syy(0,ny-1,1)+syy(1,ny-1,1))*u(0,ny-2,0)
											 +0.25*(szz(0,ny,1)+szz(1,ny,1)+szz(0,ny-1,1)+szz(1,ny-1,1))*u(0,ny-1,1)
											 +I(0,ny-1,0)*dx2)/(
											 +0.25*(sxx(1,ny,1)+sxx(1,ny-1,1)+sxx(1,ny,0)+sxx(1,ny-1,0))
											 +0.25*(syy(0,ny-1,0)+syy(1,ny-1,0)+syy(0,ny-1,1)+syy(1,ny-1,1))
											 +0.25*(szz(0,ny,1)+szz(1,ny,1)+szz(0,ny-1,1)+szz(1,ny-1,1)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		// x=0, y=0 edge
		for (k=1; k<nz-1; k++) {
			tmp = u(0,0,k);
			u(0,0,k) = (1-w)*u(0,0,k) + w*(+0.25*(sxx(1,1,k+1)+sxx(1,0,k+1)+sxx(1,1,k)+sxx(1,0,k))*u(1,0,k)
										   +0.25*(syy(0,1,k)+syy(1,1,k)+syy(0,1,k+1)+syy(1,1,k+1))*u(0,1,k)
										   +0.25*(szz(0,1,k)+szz(1,1,k)+szz(0,0,k)+szz(1,0,k))*u(0,0,k-1)
										   +0.25*(szz(0,1,k+1)+szz(1,1,k+1)+szz(0,0,k+1)+szz(1,0,k+1))*u(0,0,k+1) 
										   +I(0,0,k)*dx2)/(
										   +0.25*(sxx(1,1,k+1)+sxx(1,0,k+1)+sxx(1,1,k)+sxx(1,0,k))
										   +0.25*(syy(0,1,k)+syy(1,1,k)+syy(0,1,k+1)+syy(1,1,k+1))
										   +0.25*(szz(0,1,k)+szz(1,1,k)+szz(0,0,k)+szz(1,0,k))
										   +0.25*(szz(0,1,k+1)+szz(1,1,k+1)+szz(0,0,k+1)+szz(1,0,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		// x=0 face
		for (j=1; j<ny-1; j++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(0,j,k);
				u(0,j,k) = (1-w)*u(0,j,k) + w*(+0.25*(sxx(1,j+1,k+1)+sxx(1,j,k+1)+sxx(1,j+1,k)+sxx(1,j,k))*u(1,j,k)
											   +0.25*(syy(0,j+1,k)+syy(1,j+1,k)+syy(0,j+1,k+1)+syy(1,j+1,k+1))*u(0,j+1,k)
											   +0.25*(syy(0,j,k)+syy(1,j,k)+syy(0,j,k+1)+syy(1,j,k+1))*u(0,j-1,k)
											   +0.25*(szz(0,j+1,k)+szz(1,j+1,k)+szz(0,j,k)+szz(1,j,k))*u(0,j,k-1)
											   +0.25*(szz(0,j+1,k+1)+szz(1,j+1,k+1)+szz(0,j,k+1)+szz(1,j,k+1))*u(0,j,k+1)
											   +I(0,j,k)*dx2)/(
											   +0.25*(sxx(1,j+1,k+1)+sxx(1,j,k+1)+sxx(1,j+1,k)+sxx(1,j,k))
											   +0.25*(syy(0,j+1,k)+syy(1,j+1,k)+syy(0,j+1,k+1)+syy(1,j+1,k+1))
											   +0.25*(syy(0,j,k)+syy(1,j,k)+syy(0,j,k+1)+syy(1,j,k+1))
											   +0.25*(szz(0,j+1,k)+szz(1,j+1,k)+szz(0,j,k)+szz(1,j,k))
											   +0.25*(szz(0,j+1,k+1)+szz(1,j+1,k+1)+szz(0,j,k+1)+szz(1,j,k+1)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		// x=0, j=ny edge
		for (k=1; k<nz-1; k++) {
			tmp = u(0,ny-1,k);
			u(0,ny-1,k) = (1-w)*u(0,ny-1,k) + w*(+0.25*(sxx(1,ny-1,k+1)+sxx(1,ny-1,k+1)+sxx(1,ny,k)+sxx(1,ny-1,k))*u(1,ny-1,k)
												 +0.25*(syy(0,ny-1,k)+syy(1,ny-1,k)+syy(0,ny-1,k+1)+syy(1,ny-1,k+1))*u(0,ny-2,k)
												 +0.25*(szz(0,ny,k)+szz(1,ny,k)+szz(0,ny-1,k)+szz(1,ny-1,k))*u(0,ny-1,k-1)
												 +0.25*(szz(0,ny,k+1)+szz(1,ny,k+1)+szz(0,ny-1,k+1)+szz(1,ny-1,k+1))*u(0,ny-1,k+1)
												 +I(0,ny-1,k)*dx2)/(
												 +0.25*(sxx(1,ny-1,k+1)+sxx(1,ny-1,k+1)+sxx(1,ny,k)+sxx(1,ny-1,k))
												 +0.25*(syy(0,ny-1,k)+syy(1,ny-1,k)+syy(0,ny-1,k+1)+syy(1,ny-1,k+1))
												 +0.25*(szz(0,ny,k)+szz(1,ny,k)+szz(0,ny-1,k)+szz(1,ny-1,k))
												 +0.25*(szz(0,ny,k+1)+szz(1,ny,k+1)+szz(0,ny-1,k+1)+szz(1,ny-1,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		// x=0, y=0, z=nz corner
		tmp = u(0,0,nz-1);
		u(0,0,nz-1) = (1-w)*u(0,0,nz-1) + w*(+0.25*(sxx(1,1,nz)+sxx(1,0,nz)+sxx(1,1,nz-1)+sxx(1,0,nz-1))*u(1,0,nz-1)
											 +0.25*(syy(0,1,nz-1)+syy(1,1,nz-1)+syy(0,1,nz)+syy(1,1,nz))*u(0,1,nz-1)
											 +0.25*(szz(0,1,nz-1)+szz(1,1,nz-1)+szz(0,0,nz-1)+szz(1,0,nz-1))*u(0,0,nz-2)
											 +I(0,0,nz-1)*dx2)/(
											 +0.25*(sxx(1,1,nz)+sxx(1,0,nz)+sxx(1,1,nz-1)+sxx(1,0,nz-1))
											 +0.25*(syy(0,1,nz-1)+syy(1,1,nz-1)+syy(0,1,nz)+syy(1,1,nz))
											 +0.25*(szz(0,1,nz-1)+szz(1,1,nz-1)+szz(0,0,nz-1)+szz(1,0,nz-1)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		// x=0, z=nz edge
		for (j=1; j<ny-1; j++) {
			tmp = u(0,j,nz-1);
			u(0,j,nz-1) = (1-w)*u(0,j,nz-1) + w*(+0.25*(sxx(1,j+1,nz)+sxx(1,j,nz)+sxx(1,j+1,nz-1)+sxx(1,j,nz-1))*u(1,j,nz-1)
												 +0.25*(syy(0,j+1,nz-1)+syy(1,j+1,nz-1)+syy(0,j+1,nz)+syy(1,j+1,nz))*u(0,j+1,nz-1)
												 +0.25*(syy(0,j,nz-1)+syy(1,j,nz-1)+syy(0,j,nz)+syy(1,j,nz))*u(0,j-1,nz-1)
												 +0.25*(szz(0,j+1,nz-1)+szz(1,j+1,nz-1)+szz(0,j,nz-1)+szz(1,j,nz-1))*u(0,j,nz-2)
												 +I(0,j,nz-1)*dx2)/(
												 +0.25*(sxx(1,j+1,nz)+sxx(1,j,nz)+sxx(1,j+1,nz-1)+sxx(1,j,nz-1))
												 +0.25*(syy(0,j+1,nz-1)+syy(1,j+1,nz-1)+syy(0,j+1,nz)+syy(1,j+1,nz))
												 +0.25*(syy(0,j,nz-1)+syy(1,j,nz-1)+syy(0,j,nz)+syy(1,j,nz))
												 +0.25*(szz(0,j+1,nz-1)+szz(1,j+1,nz-1)+szz(0,j,nz-1)+szz(1,j,nz-1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		// x=0, y=ny, z=nz corner
		tmp = u(0,ny-1,nz-1);
		u(0,ny-1,nz-1) = (1-w)*u(0,ny-1,nz-1) + w*(+0.25*(sxx(1,ny,nz)+sxx(1,ny-1,nz)+sxx(1,ny,nz-1)+sxx(1,ny-1,nz-1))*u(1,ny-1,nz-1)
												   +0.25*(syy(0,ny-1,nz-1)+syy(1,ny-1,nz-1)+syy(0,ny-1,nz)+syy(1,ny,nz))*u(0,ny-2,nz-1)
												   +0.25*(szz(0,ny,nz-1)+szz(1,ny,nz-1)+szz(0,ny-1,nz-1)+szz(1,ny-1,nz-1))*u(0,ny-1,nz-2)
												   +I(0,ny-1,nz-1)*dx2)/(
												   +0.25*(sxx(1,ny,nz)+sxx(1,ny-1,nz)+sxx(1,ny,nz-1)+sxx(1,ny-1,nz-1))
												   +0.25*(syy(0,ny-1,nz-1)+syy(1,ny-1,nz-1)+syy(0,ny-1,nz)+syy(1,ny,nz))
												   +0.25*(szz(0,ny,nz-1)+szz(1,ny,nz-1)+szz(0,ny-1,nz-1)+szz(1,ny-1,nz-1)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		// j=0, z=0 edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,0,0);
			u(i,0,0) = (1-w)*u(i,0,0) + w*(+0.25*(sxx(i+1,1,1)+sxx(i+1,0,1)+sxx(i+1,1,0)+sxx(i+1,0,0))*u(i+1,0,0)
										   +0.25*(syy(i,1,0)+syy(i+1,1,0)+syy(i,1,1)+syy(i+1,1,1))*u(i,1,0)
										   +0.25*(sxx(i,1,0)+sxx(i,0,0)+sxx(i,1,1)+sxx(i,0,1))*u(i-1,0,0)
										   +0.25*(szz(i,1,1)+szz(i+1,1,1)+szz(i,0,1)+szz(i+1,0,1))*u(i,0,1)
										   +I(i,0,0)*dx2)/(
										   +0.25*(sxx(i+1,1,1)+sxx(i+1,0,1)+sxx(i+1,1,0)+sxx(i+1,0,0))
										   +0.25*(syy(i,1,0)+syy(i+1,1,0)+syy(i,1,1)+syy(i+1,1,1))
										   +0.25*(sxx(i,1,0)+sxx(i,0,0)+sxx(i,1,1)+sxx(i,0,1))
										   +0.25*(szz(i,1,1)+szz(i+1,1,1)+szz(i,0,1)+szz(i+1,0,1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}
		
		int j=0; // front face
		for (i=1; i<nx-1; i++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
											   +I(i,j,k)*dx2)/(
											   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		j=0,k=nz-1; // front top edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)										   
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		k=0; // bottom face
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
											   +I(i,j,k)*dx2)/(
											   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		j=ny-1, k=0; // back bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		// bulk
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				for (k=1; k<nz-1; k++) {
					tmp = u(i,j,k);
					u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
												   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
												   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
												   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
												   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
												   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
												   +I(i,j,k)*dx2)/(
												   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
												   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
												   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
												   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
												   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
												   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
					if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
					else {
						diff = u(i,j,k) - tmp;
						err += diff*diff;
						sum += tmp*tmp;
					}
				}
			}
		}

		k=nz-1; // top face
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
											   +I(i,j,k)*dx2)/(
											   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		j=ny-1; // back face
		for (i=1; i<nx-1; i++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
											   +I(i,j,k)*dx2)/(
											   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		j=ny-1, k=nz-1; // back top edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=nx-1, j=0, k=0; // bottom front right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
									   +I(i,j,k)*dx2)/(
									   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, k=0; // bottom right edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   +I(i,j,k)*dx2)/(
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}
		
		i=nx-1, j=ny-1, k=0; // bottom back right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
									   +I(i,j,k)*dx2)/(
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, j=0; // front right edge
		for (k=1; k<nz-1; k++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   +I(i,j,k)*dx2)/(
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=nx-1, j=0, k=nz-1; // front right top corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
									   +I(i,j,k)*dx2)/(
									   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1; // right face
		for (j=1; j<ny-1; j++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
											   +I(i,j,k)*dx2)/(
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		i=nx-1, j=ny-1; // back right edge
		for (k=1; k<nz-1; k++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=nx-1, k=nz-1; // right top edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   +I(i,j,k)*dx2)/(
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=nx-1, j=ny-1, k=nz-1; // top right back corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
									   +I(i,j,k)*dx2)/(
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		if (fabs(sum) != 0.) {
	        err = sqrt(err/sum);
	        if (err<double(tol)) break;
	    }  
	}

	if (err>double(tol)) return_val=double(err);
	else return_val = double(t);
	'''

sor_3d_inh_xy = '''
	# line 13 "sor_poisson_3d_inh_xy.py"
	double tmp, diff, err=1.0, sum, dx2 = double(_dx2);
	int t=0, i, j, k;

	for (t=0; t<int(nmax); t++) {
		err = 0.0;
		sum = 0.0;

	    i=0, j=0, k=0; // left front bottom corner neumann boundary conditions with 1st order approximation
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
									   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
									   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
									   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i,j+1,k)-u(i+1,j+1,k)-u(i+1,j,k))
									   +I(i,j,k)*dx2)/(
									   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) + 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
									   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=0, k=0; // left bottom edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   +0.0625*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))
										   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i,j+1,k)-u(i+1,j+1,k)-u(i+1,j,k))
										   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) - 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=0, j=ny-1, k=0; // left back bottom corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
									   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i+1,j,k)-u(i,j-1,k)-u(i+1,j-1,k))
									   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j,k))
									   +I(i,j,k)*dx2)/(
									   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) - 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) - 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=0, j=0; // front left edge
		for (k=1; k<nz-1; k++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
										   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i,j+1,k)-u(i+1,j+1,k)-u(i+1,j,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) + 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=0; // left face
		for (j=1; j<ny-1; j++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
											   +0.0625*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))
											   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i,j+1,k)-u(i+1,j+1,k)-u(i+1,j,k))
											   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i,j-1,k))
											   +I(i,j,k)*dx2)/(
											   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) - 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		i=0, j=ny-1; // left back edge
		for (k=1; k<nz-1; k++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i+1,j,k)-u(i,j-1,k)-u(i+1,j-1,k))
										   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i,j-1,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) - 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) - 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=0, j=0, k=nz-1; // left front top corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
									   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
									   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
									   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i,j+1,k)-u(i+1,j+1,k)-u(i+1,j,k))
									   +I(i,j,k)*dx2)/(
									   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) + 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
									   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=0, k=nz-1; // left top edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   +0.0625*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))
										   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i,j+1,k)-u(i+1,j+1,k)-u(i+1,j,k))
										   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i,j-1,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) - 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=0, j=ny-1, k=nz-1; // left back top corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
									   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i+1,j,k)-u(i,j-1,k)-u(i+1,j-1,k))
									   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i,j-1,k))
									   +I(i,j,k)*dx2)/(
									   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) - 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) - 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		j=0, k=0; // front bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
										   -0.0625*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i+1,j+1,k)-u(i+1,j,k))
										   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) + 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}
		
		int j=0; // front face
		for (i=1; i<nx-1; i++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
											   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
											   -0.0625*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i+1,j+1,k)-u(i+1,j,k))
											   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
											   +I(i,j,k)*dx2)/(
											   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) + 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		j=0,k=nz-1; // front top edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)										   
										   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
										   -0.0625*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i+1,j+1,k)-u(i+1,j,k))
										   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) + 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		k=0; // bottom face
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
											   +0.0625*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))
											   -0.0625*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i+1,j+1,k)-u(i+1,j,k))
											   +0.0625*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j+1,k)-u(i,j+1,k))
											   -0.0625*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j-1,k)-u(i-1,j,k))
											   +I(i,j,k)*dx2)/(
											   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		j=ny-1, k=0; // back bottom edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i+1,j,k))
										   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   -0.0625*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j-1,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) + 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		// bulk
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				for (k=1; k<nz-1; k++) {
					tmp = u(i,j,k);
					u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
												   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
												   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
												   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
												   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
												   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
												   +0.0625*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))
												   -0.0625*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i+1,j+1,k)-u(i+1,j,k))
												   +0.0625*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j+1,k)-u(i,j+1,k))
												   -0.0625*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j-1,k)-u(i-1,j,k))
												   +I(i,j,k)*dx2)/(
												   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
												   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
												   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
												   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
												   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
												   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
					if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
					else {
						diff = u(i,j,k) - tmp;
						err += diff*diff;
						sum += tmp*tmp;
					}
				}
			}
		}

		k=nz-1; // top face
		for (i=1; i<nx-1; i++) {
			for (j=1; j<ny-1; j++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
											   +0.0625*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))
											   -0.0625*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i+1,j+1,k)-u(i+1,j,k))
											   +0.0625*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j+1,k)-u(i,j+1,k))
											   -0.0625*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j-1,k)-u(i-1,j,k))
											   +I(i,j,k)*dx2)/(
											   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		j=ny-1; // back face
		for (i=1; i<nx-1; i++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
											   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i+1,j,k)-u(i,j-1,k)-u(i+1,j-1,k))
											   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j,k))
											   -0.0625*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j-1,k)-u(i-1,j,k))
											   +I(i,j,k)*dx2)/(
											   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) - 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		j=ny-1, k=nz-1; // back top edge
		for (i=1; i<nx-1; i++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k))*u(i+1,j,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   +0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))*(u(i+1,j,k)-u(i,j-1,k)-u(i+1,j-1,k))
										   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j,k))
										   -0.0625*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i+1,j-1,k)+u(i+1,j,k)-u(i-1,j-1,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i+1,j+1,k+1)+sxx(i+1,j,k+1)+sxx(i+1,j+1,k)+sxx(i+1,j,k)) - 0.125*(sxy(i+1,j+1,k+1)+sxy(i+1,j,k+1)+sxy(i+1,j+1,k)+sxy(i+1,j,k))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=nx-1, j=0, k=0; // bottom front right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
									   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i,j+1,k))
									   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
									   +I(i,j,k)*dx2)/(
									   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, k=0; // bottom right edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i,j+1,k))
										   +0.0625*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i,j-1,k)-u(i-1,j-1,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) + 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}
		
		i=nx-1, j=ny-1, k=0; // bottom back right corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
									   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j,k))
									   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i,j-1,k)-u(i-1,j-1,k)-u(i-1,j,k))
									   +I(i,j,k)*dx2)/(
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) + 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
									   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1, j=0; // front right edge
		for (k=1; k<nz-1; k++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i,j+1,k))
										   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=nx-1, j=0, k=nz-1; // front right top corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
									   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i,j+1,k))
									   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j,k)-u(i-1,j+1,k)-u(i,j+1,k))
									   +I(i,j,k)*dx2)/(
									   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		i=nx-1; // right face
		for (j=1; j<ny-1; j++) {
			for (k=1; k<nz-1; k++) {
				tmp = u(i,j,k);
				u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
											   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i,j+1,k))
											   +0.0625*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j+1,k)-u(i,j+1,k))
											   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i,j-1,k)-u(i-1,j-1,k)-u(i-1,j,k))
											   +I(i,j,k)*dx2)/(
											   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
											   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
											   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) + 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
											   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
											   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
				if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
				else {
					diff = u(i,j,k) - tmp;
					err += diff*diff;
					sum += tmp*tmp;
				}
			}
		}

		i=nx-1, j=ny-1; // back right edge
		for (k=1; k<nz-1; k++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1))*u(i,j,k+1)
										   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j,k))
										   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i,j-1,k)-u(i-1,j-1,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) + 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))
										   +0.25*(szz(i,j+1,k+1)+szz(i+1,j+1,k+1)+szz(i,j,k+1)+szz(i+1,j,k+1)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=nx-1, k=nz-1; // right top edge
		for (j=1; j<ny-1; j++) {
			tmp = u(i,j,k);
			u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1))*u(i,j+1,k)
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
										   -0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))*(u(i-1,j+1,k)+u(i-1,j,k)-u(i,j+1,k))
										   +0.0625*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j+1,k)-u(i,j+1,k))
										   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i,j-1,k)-u(i-1,j-1,k)-u(i-1,j,k))
										   +I(i,j,k)*dx2)/(
										   +0.25*(syy(i,j+1,k)+syy(i+1,j+1,k)+syy(i,j+1,k+1)+syy(i+1,j+1,k+1)) - 0.125*(sxy(i,j+1,k)+sxy(i+1,j+1,k)+sxy(i,j+1,k+1)+sxy(i+1,j+1,k+1))
										   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))
										   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) + 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
										   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
			if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
			else {
				diff = u(i,j,k) - tmp;
				err += diff*diff;
				sum += tmp*tmp;
			}
		}

		i=nx-1, j=ny-1, k=nz-1; // top right back corner
		tmp = u(i,j,k);
		u(i,j,k) = (1-w)*u(i,j,k) + w*(+0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1))*u(i-1,j,k)
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1))*u(i,j-1,k)
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k))*u(i,j,k-1)
									   +0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))*(u(i-1,j-1,k)+u(i,j-1,k)-u(i-1,j,k))
									   -0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))*(u(i,j-1,k)-u(i-1,j-1,k)-u(i-1,j,k))
									   +I(i,j,k)*dx2)/(
									   +0.25*(sxx(i,j+1,k)+sxx(i,j,k)+sxx(i,j+1,k+1)+sxx(i,j,k+1)) + 0.125*(sxy(i,j+1,k)+sxy(i,j,k)+sxy(i,j+1,k+1)+sxy(i,j,k+1))
									   +0.25*(syy(i,j,k)+syy(i+1,j,k)+syy(i,j,k+1)+syy(i+1,j,k+1)) + 0.125*(sxy(i,j,k)+sxy(i+1,j,k)+sxy(i,j,k+1)+sxy(i+1,j,k+1))
									   +0.25*(szz(i,j+1,k)+szz(i+1,j+1,k)+szz(i,j,k)+szz(i+1,j,k)));
		if (u(i,j,k) != u(i,j,k)) u(i,j,k) = 0;
		else {
			diff = u(i,j,k) - tmp;
			err += diff*diff;
			sum += tmp*tmp;
		}

		if (fabs(sum) != 0.) {
	        err = sqrt(err/sum);
	        if (err<double(tol)) break;
	    }  
	}

	if (err>double(tol)) return_val=double(err);
	else return_val = double(t);
	'''