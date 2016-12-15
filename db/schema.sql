/*
    Schema description for The Orchestra database.
*/

DROP TABLE IF EXISTS phase3_products;
CREATE TABLE phase3_products (
    arcfile text PRIMARY KEY,
    object text,
    ra numeric,
    dec numeric,
    wavelength text,
    snr numeric,
    resolution integer,
    instrument text,
    date_obs text,
    exptime numeric,
    program_id text,
    origfile text,
    dataset text
);

# Create Q3C indices
CREATE INDEX ON phase3_products (q3c_ang2ipix(ra, dec));
CLUSTER phase3_products_q3c_ang2ipix_idx ON phase3_products;
ANALYZE phase3_products;

/*
    Python code to generate schema for obs table:

        from astropy.io import fits

        image = fits.open("HARPS.2013-07-28T07:30:15.108_bis_G2_A.fits")

        translate_types = {int: "integer", float: "numeric", str: "text", bool: "boolean"}

        for k, v in image[0].header.items():
            # Truncate the ESO shit.
            key = k[4:] if k[:4].startswith("ESO ") else k
            key = key.replace(" ", "_").replace("-", "_").lower()

            # Ignore this garbage.
            if key.startswith("drs_cal_"): continue

            dtype = translate_types[type(v)]
            print("  {} {},".format(key, dtype))
*/

DROP TABLE IF EXISTS obs;
CREATE TABLE obs (
    date_obs text PRIMARY KEY,
    filename text,
    simple boolean,
    bitpix integer,
    naxis integer,
    naxis1 integer,
    extend boolean,
    cdelt1 numeric,
    drs_bis_rv numeric,
    drs_bis_span numeric,
    crval1 numeric,
    drs_ccf_rvc numeric not null,
    drs_ccf_contrast numeric,
    drs_ccf_fwhm numeric,
    drs_ccf_rv numeric,
    ctype1 text,
    drs_ccf_lines integer,
    drs_ccf_maxcpp integer,
    drs_ccf_mask text,
    drs_ccf_noise numeric,
    origin text,
    date text,
    telescop text,
    instrume text,
    object text,
    ra numeric,
    dec numeric,
    equinox numeric,
    radecsys text,
    exptime numeric,
    mjd_obs numeric not null,
    utc numeric,
    lst numeric,
    pi_coi text,
    ada_absrot_end numeric,
    ada_absrot_start numeric,
    ada_guid_dec numeric,
    ada_guid_ra numeric,
    ada_guid_status text,
    ada_posang numeric,
    det_bits integer,
    det_chips integer,
    det_date text,
    det_dec numeric,
    det_did text,
    det_exp_no integer,
    det_exp_rdttime numeric,
    det_exp_type text,
    det_exp_xfertim numeric,
    det_fram_id integer,
    det_fram_type text,
    det_id text,
    det_name text,
    det_outputs integer,
    det_outref integer,
    det_ra numeric,
    det_read_clock text,
    det_read_mode text,
    det_read_nfram integer,
    det_read_speed text,
    det_shut_id text,
    det_shut_tmclos numeric,
    det_shut_tmopen numeric,
    det_shut_type text,
    det_sofw_mode text,
    det_tele_int numeric,
    det_tele_no integer,
    det_tlm1_end numeric,
    det_tlm1_id text,
    det_tlm1_name text,
    det_tlm1_start numeric,
    det_tlm2_end numeric,
    det_tlm2_id text,
    det_tlm2_name text,
    det_tlm2_start numeric,
    det_tlm3_end numeric,
    det_tlm3_id text,
    det_tlm3_name text,
    det_tlm3_start numeric,
    det_win1_binx integer,
    det_win1_biny integer,
    det_win1_dit1 numeric,
    det_win1_dktm numeric,
    det_win1_ndit integer,
    det_win1_nx integer,
    det_win1_ny integer,
    det_win1_st boolean,
    det_win1_strx integer,
    det_win1_stry integer,
    det_win1_uit1 numeric,
    det_windows integer,
    dpr_catg text,
    dpr_tech text,
    dpr_type text,
    ins_adc1_dec numeric,
    ins_adc1_end numeric,
    ins_adc1_mode text,
    ins_adc1_ra numeric,
    ins_adc2_dec numeric,
    ins_adc2_end numeric,
    ins_adc2_mode text,
    ins_adc2_ra numeric,
    ins_adcs_name text,
    ins_adcs_no integer,
    ins_date text,
    ins_det1_ctmax numeric,
    ins_det1_ctmean numeric,
    ins_det1_ctmin numeric,
    ins_det1_ctrms numeric,
    ins_det1_cttot numeric,
    ins_det1_id text,
    ins_det1_name text,
    ins_det1_offdrk numeric,
    ins_det1_offsky numeric,
    ins_det1_tmmean numeric,
    ins_det1_uit numeric,
    ins_det2_ctmax numeric,
    ins_det2_ctmean numeric,
    ins_det2_ctmin numeric,
    ins_det2_ctrms numeric,
    ins_det2_cttot numeric,
    ins_det2_id text,
    ins_det2_name text,
    ins_det2_offdrk numeric,
    ins_det2_offsky numeric,
    ins_det2_tmmean numeric,
    ins_det2_uit numeric,
    ins_did text,
    ins_hefs_st boolean,
    ins_id text,
    ins_lamp6_id text,
    ins_lamp6_name text,
    ins_lamp6_st boolean,
    ins_lamp7_st boolean,
    ins_mirr1_id text,
    ins_mirr1_name text,
    ins_mirr1_no integer,
    ins_mode text,
    ins_opti2_id text,
    ins_opti2_name text,
    ins_opti2_no integer,
    ins_opti2_type text,
    ins_opti3_id text,
    ins_opti3_name text,
    ins_opti3_no integer,
    ins_opti3_type text,
    ins_opti4_id text,
    ins_opti4_name text,
    ins_opti4_no integer,
    ins_opti4_swsim boolean,
    ins_opti4_type text,
    ins_opti5_id text,
    ins_opti5_name text,
    ins_opti5_no integer,
    ins_opti5_type text,
    ins_opti6_id text,
    ins_opti6_name text,
    ins_opti6_no integer,
    ins_opti6_type text,
    ins_opti7_id text,
    ins_opti7_name text,
    ins_opti7_no integer,
    ins_opti7_type text,
    ins_path text,
    ins_rot1_enc integer,
    ins_rot1_pos numeric,
    ins_sens1_detcoef numeric,
    ins_sens1_grad numeric,
    ins_sens1_id text,
    ins_sens1_lrconst numeric,
    ins_sens1_lrrms numeric,
    ins_sens1_max numeric,
    ins_sens1_mean numeric,
    ins_sens1_min numeric,
    ins_sens1_name text,
    ins_sens1_rms numeric,
    ins_sens1_val numeric,
    ins_sens2_detcoef numeric,
    ins_sens2_grad numeric,
    ins_sens2_id text,
    ins_sens2_lrconst numeric,
    ins_sens2_lrrms numeric,
    ins_sens2_max numeric,
    ins_sens2_mean numeric,
    ins_sens2_min numeric,
    ins_sens2_name text,
    ins_sens2_rms numeric,
    ins_sens2_val numeric,
    ins_sens6_detcoef numeric,
    ins_sens6_grad numeric,
    ins_sens6_id text,
    ins_sens6_lrconst numeric,
    ins_sens6_lrrms numeric,
    ins_sens6_max numeric,
    ins_sens6_mean numeric,
    ins_sens6_min numeric,
    ins_sens6_name text,
    ins_sens6_rms numeric,
    ins_sens6_val numeric,
    ins_sensor2_swsim boolean,
    ins_sensor4_swsim boolean,
    ins_sensor5_swsim boolean,
    ins_sensor6_swsim boolean,
    ins_sensor7_swsim boolean,
    ins_sensor8_swsim boolean,
    ins_sensor9_swsim boolean,
    ins_shut1_id text,
    ins_shut1_name text,
    ins_shut1_st boolean,
    ins_swsim text,
    ins_temp1_detcoef numeric,
    ins_temp1_grad numeric,
    ins_temp1_id text,
    ins_temp1_lrconst numeric,
    ins_temp1_lrrms numeric,
    ins_temp1_max numeric,
    ins_temp1_mean numeric,
    ins_temp1_min numeric,
    ins_temp1_name text,
    ins_temp1_rms numeric,
    ins_temp1_val numeric,
    ins_temp10_detcoef numeric,
    ins_temp10_grad numeric,
    ins_temp10_id text,
    ins_temp10_lrconst numeric,
    ins_temp10_lrrms numeric,
    ins_temp10_max numeric,
    ins_temp10_mean numeric,
    ins_temp10_min numeric,
    ins_temp10_name text,
    ins_temp10_rms numeric,
    ins_temp10_val numeric,
    ins_temp11_detcoef numeric,
    ins_temp11_grad numeric,
    ins_temp11_id text,
    ins_temp11_lrconst numeric,
    ins_temp11_lrrms numeric,
    ins_temp11_max numeric,
    ins_temp11_mean numeric,
    ins_temp11_min numeric,
    ins_temp11_name text,
    ins_temp11_rms numeric,
    ins_temp11_val numeric,
    ins_temp12_detcoef numeric,
    ins_temp12_grad numeric,
    ins_temp12_id text,
    ins_temp12_lrconst numeric,
    ins_temp12_lrrms numeric,
    ins_temp12_max numeric,
    ins_temp12_mean numeric,
    ins_temp12_min numeric,
    ins_temp12_name text,
    ins_temp12_rms numeric,
    ins_temp12_val numeric,
    ins_temp13_detcoef numeric,
    ins_temp13_grad numeric,
    ins_temp13_id text,
    ins_temp13_lrconst numeric,
    ins_temp13_lrrms numeric,
    ins_temp13_max numeric,
    ins_temp13_mean numeric,
    ins_temp13_min numeric,
    ins_temp13_name text,
    ins_temp13_rms numeric,
    ins_temp13_val numeric,
    ins_temp14_detcoef numeric,
    ins_temp14_grad numeric,
    ins_temp14_id text,
    ins_temp14_lrconst numeric,
    ins_temp14_lrrms numeric,
    ins_temp14_max numeric,
    ins_temp14_mean numeric,
    ins_temp14_min numeric,
    ins_temp14_name text,
    ins_temp14_rms numeric,
    ins_temp14_val numeric,
    ins_temp15_detcoef numeric,
    ins_temp15_grad numeric,
    ins_temp15_id text,
    ins_temp15_lrconst numeric,
    ins_temp15_lrrms numeric,
    ins_temp15_max numeric,
    ins_temp15_mean numeric,
    ins_temp15_min numeric,
    ins_temp15_name text,
    ins_temp15_rms numeric,
    ins_temp15_val numeric,
    ins_temp16_detcoef numeric,
    ins_temp16_grad numeric,
    ins_temp16_id text,
    ins_temp16_lrconst numeric,
    ins_temp16_lrrms numeric,
    ins_temp16_max numeric,
    ins_temp16_mean numeric,
    ins_temp16_min numeric,
    ins_temp16_name text,
    ins_temp16_rms numeric,
    ins_temp16_val numeric,
    ins_temp17_detcoef numeric,
    ins_temp17_grad numeric,
    ins_temp17_id text,
    ins_temp17_lrconst numeric,
    ins_temp17_lrrms numeric,
    ins_temp17_max numeric,
    ins_temp17_mean numeric,
    ins_temp17_min numeric,
    ins_temp17_name text,
    ins_temp17_rms numeric,
    ins_temp17_val numeric,
    ins_temp18_detcoef numeric,
    ins_temp18_grad numeric,
    ins_temp18_id text,
    ins_temp18_lrconst numeric,
    ins_temp18_lrrms numeric,
    ins_temp18_max numeric,
    ins_temp18_mean numeric,
    ins_temp18_min numeric,
    ins_temp18_name text,
    ins_temp18_rms numeric,
    ins_temp18_val numeric,
    ins_temp19_id text,
    ins_temp19_name text,
    ins_temp19_val numeric,
    ins_temp2_id text,
    ins_temp2_name text,
    ins_temp2_val numeric,
    ins_temp20_id text,
    ins_temp20_name text,
    ins_temp20_val numeric,
    ins_temp21_detcoef numeric,
    ins_temp21_grad numeric,
    ins_temp21_id text,
    ins_temp21_lrconst numeric,
    ins_temp21_lrrms numeric,
    ins_temp21_max numeric,
    ins_temp21_mean numeric,
    ins_temp21_min numeric,
    ins_temp21_name text,
    ins_temp21_rms numeric,
    ins_temp21_val numeric,
    ins_temp22_detcoef numeric,
    ins_temp22_grad numeric,
    ins_temp22_id text,
    ins_temp22_lrconst numeric,
    ins_temp22_lrrms numeric,
    ins_temp22_max numeric,
    ins_temp22_mean numeric,
    ins_temp22_min numeric,
    ins_temp22_name text,
    ins_temp22_rms numeric,
    ins_temp22_val numeric,
    ins_temp23_detcoef numeric,
    ins_temp23_grad numeric,
    ins_temp23_id text,
    ins_temp23_lrconst numeric,
    ins_temp23_lrrms numeric,
    ins_temp23_max numeric,
    ins_temp23_mean numeric,
    ins_temp23_min numeric,
    ins_temp23_name text,
    ins_temp23_rms numeric,
    ins_temp23_val numeric,
    ins_temp24_detcoef numeric,
    ins_temp24_grad numeric,
    ins_temp24_id text,
    ins_temp24_lrconst numeric,
    ins_temp24_lrrms numeric,
    ins_temp24_max numeric,
    ins_temp24_mean numeric,
    ins_temp24_min numeric,
    ins_temp24_name text,
    ins_temp24_rms numeric,
    ins_temp24_val numeric,
    ins_temp25_detcoef numeric,
    ins_temp25_grad numeric,
    ins_temp25_id text,
    ins_temp25_lrconst numeric,
    ins_temp25_lrrms numeric,
    ins_temp25_max numeric,
    ins_temp25_mean numeric,
    ins_temp25_min numeric,
    ins_temp25_name text,
    ins_temp25_rms numeric,
    ins_temp25_val numeric,
    ins_temp26_detcoef numeric,
    ins_temp26_grad numeric,
    ins_temp26_id text,
    ins_temp26_lrconst numeric,
    ins_temp26_lrrms numeric,
    ins_temp26_max numeric,
    ins_temp26_mean numeric,
    ins_temp26_min numeric,
    ins_temp26_name text,
    ins_temp26_rms numeric,
    ins_temp26_val numeric,
    ins_temp27_detcoef numeric,
    ins_temp27_grad numeric,
    ins_temp27_id text,
    ins_temp27_lrconst numeric,
    ins_temp27_lrrms numeric,
    ins_temp27_max numeric,
    ins_temp27_mean numeric,
    ins_temp27_min numeric,
    ins_temp27_name text,
    ins_temp27_rms numeric,
    ins_temp27_val numeric,
    ins_temp28_detcoef numeric,
    ins_temp28_grad numeric,
    ins_temp28_id text,
    ins_temp28_lrconst numeric,
    ins_temp28_lrrms numeric,
    ins_temp28_max numeric,
    ins_temp28_mean numeric,
    ins_temp28_min numeric,
    ins_temp28_name text,
    ins_temp28_rms numeric,
    ins_temp28_val numeric,
    ins_temp29_detcoef numeric,
    ins_temp29_grad numeric,
    ins_temp29_id text,
    ins_temp29_lrconst numeric,
    ins_temp29_lrrms numeric,
    ins_temp29_max numeric,
    ins_temp29_mean numeric,
    ins_temp29_min numeric,
    ins_temp29_name text,
    ins_temp29_rms numeric,
    ins_temp29_val numeric,
    ins_temp3_detcoef numeric,
    ins_temp3_grad numeric,
    ins_temp3_id text,
    ins_temp3_lrconst numeric,
    ins_temp3_lrrms numeric,
    ins_temp3_max numeric,
    ins_temp3_mean numeric,
    ins_temp3_min numeric,
    ins_temp3_name text,
    ins_temp3_rms numeric,
    ins_temp3_val numeric,
    ins_temp30_detcoef numeric,
    ins_temp30_grad numeric,
    ins_temp30_id text,
    ins_temp30_lrconst numeric,
    ins_temp30_lrrms numeric,
    ins_temp30_max numeric,
    ins_temp30_mean numeric,
    ins_temp30_min numeric,
    ins_temp30_name text,
    ins_temp30_rms numeric,
    ins_temp30_val numeric,
    ins_temp31_detcoef numeric,
    ins_temp31_grad numeric,
    ins_temp31_id text,
    ins_temp31_lrconst numeric,
    ins_temp31_lrrms numeric,
    ins_temp31_max numeric,
    ins_temp31_mean numeric,
    ins_temp31_min numeric,
    ins_temp31_name text,
    ins_temp31_rms numeric,
    ins_temp31_val numeric,
    ins_temp32_detcoef numeric,
    ins_temp32_grad numeric,
    ins_temp32_id text,
    ins_temp32_lrconst numeric,
    ins_temp32_lrrms numeric,
    ins_temp32_max numeric,
    ins_temp32_mean numeric,
    ins_temp32_min numeric,
    ins_temp32_name text,
    ins_temp32_rms numeric,
    ins_temp32_val numeric,
    ins_temp33_detcoef numeric,
    ins_temp33_grad numeric,
    ins_temp33_id text,
    ins_temp33_lrconst numeric,
    ins_temp33_lrrms numeric,
    ins_temp33_max numeric,
    ins_temp33_mean numeric,
    ins_temp33_min numeric,
    ins_temp33_name text,
    ins_temp33_rms numeric,
    ins_temp33_val numeric,
    ins_temp34_detcoef numeric,
    ins_temp34_grad numeric,
    ins_temp34_id text,
    ins_temp34_lrconst numeric,
    ins_temp34_lrrms numeric,
    ins_temp34_max numeric,
    ins_temp34_mean numeric,
    ins_temp34_min numeric,
    ins_temp34_name text,
    ins_temp34_rms numeric,
    ins_temp34_val numeric,
    ins_temp35_detcoef numeric,
    ins_temp35_grad numeric,
    ins_temp35_id text,
    ins_temp35_lrconst numeric,
    ins_temp35_lrrms numeric,
    ins_temp35_max numeric,
    ins_temp35_mean numeric,
    ins_temp35_min numeric,
    ins_temp35_name text,
    ins_temp35_rms numeric,
    ins_temp35_val numeric,
    ins_temp36_detcoef numeric,
    ins_temp36_grad numeric,
    ins_temp36_id text,
    ins_temp36_lrconst numeric,
    ins_temp36_lrrms numeric,
    ins_temp36_max numeric,
    ins_temp36_mean numeric,
    ins_temp36_min numeric,
    ins_temp36_name text,
    ins_temp36_rms numeric,
    ins_temp36_val numeric,
    ins_temp37_detcoef numeric,
    ins_temp37_grad numeric,
    ins_temp37_id text,
    ins_temp37_lrconst numeric,
    ins_temp37_lrrms numeric,
    ins_temp37_max numeric,
    ins_temp37_mean numeric,
    ins_temp37_min numeric,
    ins_temp37_name text,
    ins_temp37_rms numeric,
    ins_temp37_val numeric,
    ins_temp38_id text,
    ins_temp38_name text,
    ins_temp38_val numeric,
    ins_temp4_detcoef numeric,
    ins_temp4_grad numeric,
    ins_temp4_id text,
    ins_temp4_lrconst numeric,
    ins_temp4_lrrms numeric,
    ins_temp4_max numeric,
    ins_temp4_mean numeric,
    ins_temp4_min numeric,
    ins_temp4_name text,
    ins_temp4_rms numeric,
    ins_temp4_val numeric,
    ins_temp40_detcoef numeric,
    ins_temp40_grad numeric,
    ins_temp40_id text,
    ins_temp40_lrconst numeric,
    ins_temp40_lrrms numeric,
    ins_temp40_max numeric,
    ins_temp40_mean numeric,
    ins_temp40_min numeric,
    ins_temp40_name text,
    ins_temp40_rms numeric,
    ins_temp40_val numeric,
    ins_temp41_detcoef numeric,
    ins_temp41_grad numeric,
    ins_temp41_id text,
    ins_temp41_lrconst numeric,
    ins_temp41_lrrms numeric,
    ins_temp41_max numeric,
    ins_temp41_mean numeric,
    ins_temp41_min numeric,
    ins_temp41_name text,
    ins_temp41_rms numeric,
    ins_temp41_val numeric,
    ins_temp42_detcoef numeric,
    ins_temp42_grad numeric,
    ins_temp42_id text,
    ins_temp42_lrconst numeric,
    ins_temp42_lrrms numeric,
    ins_temp42_max numeric,
    ins_temp42_mean numeric,
    ins_temp42_min numeric,
    ins_temp42_name text,
    ins_temp42_rms numeric,
    ins_temp42_val numeric,
    ins_temp43_detcoef numeric,
    ins_temp43_grad numeric,
    ins_temp43_id text,
    ins_temp43_lrconst numeric,
    ins_temp43_lrrms numeric,
    ins_temp43_max numeric,
    ins_temp43_mean numeric,
    ins_temp43_min numeric,
    ins_temp43_name text,
    ins_temp43_rms numeric,
    ins_temp43_val numeric,
    ins_temp44_detcoef numeric,
    ins_temp44_grad numeric,
    ins_temp44_id text,
    ins_temp44_lrconst numeric,
    ins_temp44_lrrms numeric,
    ins_temp44_max numeric,
    ins_temp44_mean numeric,
    ins_temp44_min numeric,
    ins_temp44_name text,
    ins_temp44_rms numeric,
    ins_temp44_val numeric,
    ins_temp5_detcoef numeric,
    ins_temp5_grad numeric,
    ins_temp5_id text,
    ins_temp5_lrconst numeric,
    ins_temp5_lrrms numeric,
    ins_temp5_max numeric,
    ins_temp5_mean numeric,
    ins_temp5_min numeric,
    ins_temp5_name text,
    ins_temp5_rms numeric,
    ins_temp5_val numeric,
    ins_temp6_detcoef numeric,
    ins_temp6_grad numeric,
    ins_temp6_id text,
    ins_temp6_lrconst numeric,
    ins_temp6_lrrms numeric,
    ins_temp6_max numeric,
    ins_temp6_mean numeric,
    ins_temp6_min numeric,
    ins_temp6_name text,
    ins_temp6_rms numeric,
    ins_temp6_val numeric,
    ins_temp7_detcoef numeric,
    ins_temp7_grad numeric,
    ins_temp7_id text,
    ins_temp7_lrconst numeric,
    ins_temp7_lrrms numeric,
    ins_temp7_max numeric,
    ins_temp7_mean numeric,
    ins_temp7_min numeric,
    ins_temp7_name text,
    ins_temp7_rms numeric,
    ins_temp7_val numeric,
    ins_temp8_detcoef numeric,
    ins_temp8_grad numeric,
    ins_temp8_id text,
    ins_temp8_lrconst numeric,
    ins_temp8_lrrms numeric,
    ins_temp8_max numeric,
    ins_temp8_mean numeric,
    ins_temp8_min numeric,
    ins_temp8_name text,
    ins_temp8_rms numeric,
    ins_temp8_val numeric,
    ins_temp9_detcoef numeric,
    ins_temp9_grad numeric,
    ins_temp9_id text,
    ins_temp9_lrconst numeric,
    ins_temp9_lrrms numeric,
    ins_temp9_max numeric,
    ins_temp9_mean numeric,
    ins_temp9_min numeric,
    ins_temp9_name text,
    ins_temp9_rms numeric,
    ins_temp9_val numeric,
    obs_did text,
    obs_grp text,
    obs_id integer,
    obs_name text,
    obs_pi_coi_id integer,
    obs_pi_coi_name text,
    obs_prog_id text,
    obs_start text,
    obs_targ_name text,
    obs_tplno integer,
    ocs_det1_imgname text,
    ocs_det1_naming text,
    tel_airm_end numeric,
    tel_airm_start numeric,
    tel_alt numeric,
    tel_ambi_fwhm_end numeric,
    tel_ambi_fwhm_start numeric,
    tel_ambi_pres_end numeric,
    tel_ambi_pres_start numeric,
    tel_ambi_rhum numeric,
    tel_ambi_temp numeric,
    tel_ambi_winddir numeric,
    tel_ambi_windsp numeric,
    tel_az numeric,
    tel_chop_st boolean,
    tel_date text,
    tel_did text,
    tel_dome_status text,
    tel_focu_id text,
    tel_focu_len numeric,
    tel_focu_scale numeric,
    tel_focu_value numeric,
    tel_geoelev numeric,
    tel_geolat numeric,
    tel_geolon numeric,
    tel_id text,
    tel_moon_dec numeric,
    tel_moon_ra numeric,
    tel_oper text,
    tel_parang_end numeric,
    tel_parang_start numeric,
    tel_targ_alpha numeric,
    tel_targ_coordtype text,
    tel_targ_delta numeric,
    tel_targ_epoch numeric,
    tel_targ_epochsystem text,
    tel_targ_equinox numeric,
    tel_targ_parallax numeric,
    tel_targ_pma numeric,
    tel_targ_pmd numeric,
    tel_targ_radvel numeric,
    tel_th_m1_temp numeric,
    tel_trak_status text,
    tpl_did text,
    tpl_expno integer,
    tpl_id text,
    tpl_name text,
    tpl_nexp integer,
    tpl_preseq text,
    tpl_start text,
    tpl_version text,
    origfile text,
    datasum text,
    drs_version text,
    drs_ccd_sigdet numeric,
    drs_ccd_conad numeric,
    drs_cal_th_file text,
    drs_cal_flat_file text,
    drs_blaze_file text,
    drs_spe_ext_opt integer,
    drs_spe_ext_window integer,
    drs_spe_ext_cosm numeric,
    drs_spe_ext_sn0 numeric,
    drs_spe_ext_sn1 numeric,
    drs_spe_ext_sn2 numeric,
    drs_spe_ext_sn3 numeric,
    drs_spe_ext_sn4 numeric,
    drs_spe_ext_sn5 numeric,
    drs_spe_ext_sn6 numeric,
    drs_spe_ext_sn7 numeric,
    drs_spe_ext_sn8 numeric,
    drs_spe_ext_sn9 numeric,
    drs_spe_ext_sn10 numeric,
    drs_spe_ext_sn11 numeric,
    drs_spe_ext_sn12 numeric,
    drs_spe_ext_sn13 numeric,
    drs_spe_ext_sn14 numeric,
    drs_spe_ext_sn15 numeric,
    drs_spe_ext_sn16 numeric,
    drs_spe_ext_sn17 numeric,
    drs_spe_ext_sn18 numeric,
    drs_spe_ext_sn19 numeric,
    drs_spe_ext_sn20 numeric,
    drs_spe_ext_sn21 numeric,
    drs_spe_ext_sn22 numeric,
    drs_spe_ext_sn23 numeric,
    drs_spe_ext_sn24 numeric,
    drs_spe_ext_sn25 numeric,
    drs_spe_ext_sn26 numeric,
    drs_spe_ext_sn27 numeric,
    drs_spe_ext_sn28 numeric,
    drs_spe_ext_sn29 numeric,
    drs_spe_ext_sn30 numeric,
    drs_spe_ext_sn31 numeric,
    drs_spe_ext_sn32 numeric,
    drs_spe_ext_sn33 numeric,
    drs_spe_ext_sn34 numeric,
    drs_spe_ext_sn35 numeric,
    drs_spe_ext_sn36 numeric,
    drs_spe_ext_sn37 numeric,
    drs_spe_ext_sn38 numeric,
    drs_spe_ext_sn39 numeric,
    drs_spe_ext_sn40 numeric,
    drs_spe_ext_sn41 numeric,
    drs_spe_ext_sn42 numeric,
    drs_spe_ext_sn43 numeric,
    drs_spe_ext_sn44 numeric,
    drs_spe_ext_sn45 numeric,
    drs_spe_ext_sn46 numeric,
    drs_spe_ext_sn47 numeric,
    drs_spe_ext_sn48 numeric,
    drs_spe_ext_sn49 numeric,
    drs_spe_ext_sn50 numeric,
    drs_spe_ext_sn51 numeric,
    drs_spe_ext_sn52 numeric,
    drs_spe_ext_sn53 numeric,
    drs_spe_ext_sn54 numeric,
    drs_spe_ext_sn55 numeric,
    drs_spe_ext_sn56 numeric,
    drs_spe_ext_sn57 numeric,
    drs_spe_ext_sn58 numeric,
    drs_spe_ext_sn59 numeric,
    drs_spe_ext_sn60 numeric,
    drs_spe_ext_sn61 numeric,
    drs_spe_ext_sn62 numeric,
    drs_spe_ext_sn63 numeric,
    drs_spe_ext_sn64 numeric,
    drs_spe_ext_sn65 numeric,
    drs_spe_ext_sn66 numeric,
    drs_spe_ext_sn67 numeric,
    drs_spe_ext_sn68 numeric,
    drs_spe_ext_sn69 numeric,
    drs_spe_ext_sn70 numeric,
    drs_spe_ext_sn71 numeric,
    drs_spe_nbcos0 integer,
    drs_spe_nbcos1 integer,
    drs_spe_nbcos2 integer,
    drs_spe_nbcos3 integer,
    drs_spe_nbcos4 integer,
    drs_spe_nbcos5 integer,
    drs_spe_nbcos6 integer,
    drs_spe_nbcos7 integer,
    drs_spe_nbcos8 integer,
    drs_spe_nbcos9 integer,
    drs_spe_nbcos10 integer,
    drs_spe_nbcos11 integer,
    drs_spe_nbcos12 integer,
    drs_spe_nbcos13 integer,
    drs_spe_nbcos14 integer,
    drs_spe_nbcos15 integer,
    drs_spe_nbcos16 integer,
    drs_spe_nbcos17 integer,
    drs_spe_nbcos18 integer,
    drs_spe_nbcos19 integer,
    drs_spe_nbcos20 integer,
    drs_spe_nbcos21 integer,
    drs_spe_nbcos22 integer,
    drs_spe_nbcos23 integer,
    drs_spe_nbcos24 integer,
    drs_spe_nbcos25 integer,
    drs_spe_nbcos26 integer,
    drs_spe_nbcos27 integer,
    drs_spe_nbcos28 integer,
    drs_spe_nbcos29 integer,
    drs_spe_nbcos30 integer,
    drs_spe_nbcos31 integer,
    drs_spe_nbcos32 integer,
    drs_spe_nbcos33 integer,
    drs_spe_nbcos34 integer,
    drs_spe_nbcos35 integer,
    drs_spe_nbcos36 integer,
    drs_spe_nbcos37 integer,
    drs_spe_nbcos38 integer,
    drs_spe_nbcos39 integer,
    drs_spe_nbcos40 integer,
    drs_spe_nbcos41 integer,
    drs_spe_nbcos42 integer,
    drs_spe_nbcos43 integer,
    drs_spe_nbcos44 integer,
    drs_spe_nbcos45 integer,
    drs_spe_nbcos46 integer,
    drs_spe_nbcos47 integer,
    drs_spe_nbcos48 integer,
    drs_spe_nbcos49 integer,
    drs_spe_nbcos50 integer,
    drs_spe_nbcos51 integer,
    drs_spe_nbcos52 integer,
    drs_spe_nbcos53 integer,
    drs_spe_nbcos54 integer,
    drs_spe_nbcos55 integer,
    drs_spe_nbcos56 integer,
    drs_spe_nbcos57 integer,
    drs_spe_nbcos58 integer,
    drs_spe_nbcos59 integer,
    drs_spe_nbcos60 integer,
    drs_spe_nbcos61 integer,
    drs_spe_nbcos62 integer,
    drs_spe_nbcos63 integer,
    drs_spe_nbcos64 integer,
    drs_spe_nbcos65 integer,
    drs_spe_nbcos66 integer,
    drs_spe_nbcos67 integer,
    drs_spe_nbcos68 integer,
    drs_spe_nbcos69 integer,
    drs_spe_nbcos70 integer,
    drs_spe_nbcos71 integer,
    drs_berv numeric,
    drs_bjd numeric,
    drs_bervmx numeric,
    drs_drift_ref_spe text,
    drs_drift_spe_rv numeric,
    drs_drift_nbcos integer,
    drs_drift_rflux numeric,
    drs_drift_nbordkill integer,
    drs_drift_noise numeric,
    drs_drift_ref_ccf text,
    drs_drift_ref_rv numeric,
    drs_drift_ccf_rv numeric,
    drs_drift_rv_used numeric,
    drs_drift_algo text,
    drs_drift_qc text,
    drs_dvrms numeric not null,
    drs_flux_corr_min numeric,
    drs_flux_corr_max numeric,
    drs_flux_corr_coeff0 numeric,
    drs_flux_corr_coeff1 numeric,
    drs_flux_corr_coeff2 numeric,
    drs_flux_corr_coeff3 numeric,
    drs_flux_corr_coeff4 numeric,
    drs_flux_corr_coeff5 numeric
);


CREATE INDEX ON obs (q3c_ang2ipix(ra, dec));
CLUSTER obs_q3c_ang2ipix_idx ON obs;
ANALYZE obs;


DROP TABLE IF EXISTS stellar_activity;
CREATE TABLE stellar_activity (
    date_obs text PRIMARY KEY,
    filename text not null,
    git_hash text not null,
    s_hk numeric,
    e_s_hk numeric
);
