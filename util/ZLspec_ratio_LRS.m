function Iratio = ZLspec_ratio_LRS(band,base_band)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the ZL intensity at 'band' relative to base_band: I_band/I_base_band
%
% Matsuura et al. 2017 Fig10 data (Sep 5th, 2017 email)
% wavelength(um), relative intensity normalized at 1.25um, error
%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=[
0.83		1.268	0.0189
0.86		1.239	0.0164
0.90		1.204	0.0138
0.95		1.169	0.0116
1.00		1.204	0.0091
1.05		1.165	0.0091
1.11		1.092	0.0077
1.18		1.050	0.0048
1.25		1.000	0.0033
1.33		0.931	0.0067
1.42		0.857	0.0088
1.51		0.829	0.0062
1.60		0.764	0.0087
1.70		0.651	0.0069
];

wl_arr = data(:,1);
I_arr = data(:,2);

Ibase = spline(wl_arr,I_arr,base_band);
Iband = spline(wl_arr,I_arr,band);

Iratio = Iband./Ibase;
end