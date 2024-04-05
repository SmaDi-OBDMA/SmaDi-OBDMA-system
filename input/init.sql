--
-- PostgreSQL database dump
--

-- Dumped from database version 14.11 (Ubuntu 14.11-0ubuntu0.22.04.1)
-- Dumped by pg_dump version 14.11 (Ubuntu 14.11-0ubuntu0.22.04.1)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: plpython3u; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS plpython3u WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpython3u; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpython3u IS 'PL/Python3U untrusted procedural language';


--
-- Name: dblink; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS dblink WITH SCHEMA public;


--
-- Name: EXTENSION dblink; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION dblink IS 'connect to other postgresQL databases from within a database';


--
-- Name: uuid-ossp; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS "uuid-ossp" WITH SCHEMA public;


--
-- Name: EXTENSION "uuid-ossp"; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION "uuid-ossp" IS 'generate universally unique identifiers (UUIDs)';


--
-- Name: aggregat_type_creating_intern; Type: TYPE; Schema: public; Owner: mena
--

CREATE TYPE public.aggregat_type_creating_intern AS (
	concr_param_id_pre character varying,
	mod_meas_id_pre character varying,
	param_ids character varying[],
	mat_sample_id character varying,
	value double precision[],
	mod_id character varying
);


ALTER TYPE public.aggregat_type_creating_intern OWNER TO mena;

--
-- Name: interpol_type; Type: TYPE; Schema: public; Owner: mena
--

CREATE TYPE public.interpol_type AS (
	x real,
	y real
);


ALTER TYPE public.interpol_type OWNER TO mena;

--
-- Name: matrix_filter_agg_type; Type: TYPE; Schema: public; Owner: mena
--

CREATE TYPE public.matrix_filter_agg_type AS (
	step integer,
	matr character varying
);


ALTER TYPE public.matrix_filter_agg_type OWNER TO mena;

--
-- Name: operator_types; Type: TYPE; Schema: public; Owner: mena
--

CREATE TYPE public.operator_types AS ENUM (
    '=',
    '!=',
    '<',
    '>',
    '>=',
    '<='
);


ALTER TYPE public.operator_types OWNER TO mena;

--
-- Name: param_meas_time_type_creating_intern; Type: TYPE; Schema: public; Owner: mena
--

CREATE TYPE public.param_meas_time_type_creating_intern AS (
	concr_param_id_pre character varying,
	mod_meas_id_pre character varying,
	param_ids character varying[],
	mat_sample_id character varying,
	meas_time integer,
	value double precision[],
	mod_id character varying
);


ALTER TYPE public.param_meas_time_type_creating_intern OWNER TO mena;

--
-- Name: param_values_type_input; Type: TYPE; Schema: public; Owner: mena
--

CREATE TYPE public.param_values_type_input AS (
	param_id character varying[],
	ar character varying[],
	mat_sample_id character varying[],
	mod_id character varying,
	meas_time integer,
	mod_type character varying
);


ALTER TYPE public.param_values_type_input OWNER TO mena;

--
-- Name: param_values_type_output; Type: TYPE; Schema: public; Owner: mena
--

CREATE TYPE public.param_values_type_output AS (
	concr_param_id character varying,
	mod_meas_id character varying,
	param_id character varying,
	mat_sample_id character varying,
	meas_time integer,
	value double precision,
	mod_id character varying
);


ALTER TYPE public.param_values_type_output OWNER TO mena;

--
-- Name: param_values_type_output_matrix; Type: TYPE; Schema: public; Owner: mena
--

CREATE TYPE public.param_values_type_output_matrix AS (
	concr_param_id character varying,
	mod_meas_id character varying,
	param_id character varying,
	mat_sample_id character varying,
	meas_time integer,
	value double precision[],
	mod_id character varying
);


ALTER TYPE public.param_values_type_output_matrix OWNER TO mena;

--
-- Name: aggregat_simple(character varying, double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.aggregat_simple(_function character varying, _params double precision[]) RETURNS double precision
    LANGUAGE plpgsql STABLE
    AS $$
declare  
_out double precision;
_sql text := 'SELECT cast(%I(x) as double precision) as answer from unnest(cast(%L as double precision[])) as x';
Begin   

EXECUTE format(_sql, _function,_params) INTO _out;
RETURN _out;
END
$$;


ALTER FUNCTION public.aggregat_simple(_function character varying, _params double precision[]) OWNER TO mena;

--
-- Name: array_contained(character varying[], character varying[], integer); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.array_contained(array1 character varying[], array2 character varying[], number_same_el integer) RETURNS boolean
    LANGUAGE plpgsql STABLE
    AS $$                                              
 declare                                          
 count int = 0;
      el varchar;
 Begin 
      FOR el IN
      SELECT * FROM unnest(array2)
      LOOP
            IF el =any(array1) THEN
                  count := count +1;
            END IF;
      END LOOP;
    
      RETURN count >= number_same_el;
 END                                             
 $$;


ALTER FUNCTION public.array_contained(array1 character varying[], array2 character varying[], number_same_el integer) OWNER TO mena;

--
-- Name: calc_flux_density(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_flux_density(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql
    AS $$
declare 
Begin
	RETURN ARRAY[_params[1]*_params[2]*_params[3]];
END$$;


ALTER FUNCTION public.calc_flux_density(_params double precision[]) OWNER TO mena;

--
-- Name: calc_mag_energy_density(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_mag_energy_density(ar double precision[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$
import sympy as sy

matrix = ar
return_array = []

for mrow in range(len(matrix)):
    if mrow == 0:
        return_array.append([matrix[mrow][0], matrix[mrow][0] ,0])
    else:
        if (matrix[mrow-1][3]==matrix[mrow][3]):
            oldvalue = return_array[mrow-1][2]
            diffH = matrix[mrow][1] - matrix[mrow-1][1]
            diffJ = matrix[mrow][2] + matrix[mrow-1][2]
            diff = diffH*diffJ/2
           
            return_array.append([matrix[mrow][0], matrix[mrow][0], oldvalue+diff])
        else:
            return_array.append([matrix[mrow][0], matrix[mrow][0], 0])

return return_array
$$;


ALTER FUNCTION public.calc_mag_energy_density(ar double precision[]) OWNER TO mena;

--
-- Name: calc_mag_polarization(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_mag_polarization(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql
    AS $$
declare  
Begin   
	RETURN ARRAY[_params[1]*_params[2]];
END$$;


ALTER FUNCTION public.calc_mag_polarization(_params double precision[]) OWNER TO mena;

--
-- Name: calc_magn_stress(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_magn_stress(ar double precision[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$
import sympy as sy


axis_one = ar[0]
axis_two = ar[1]

if (int(axis_one[2]) == 1 and int(axis_two[2]) == 0):
    calc = (axis_two[1]-axis_one[1])/axis_one[3]
elif (int(axis_one[2]) == 0 and int(axis_two[2]) == 1):
    calc = (axis_one[1]-axis_two[1])/axis_one[3]
else:
    return 

return [axis_one[0], calc]

$$;


ALTER FUNCTION public.calc_magn_stress(ar double precision[]) OWNER TO mena;

--
-- Name: calc_magnetization(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_magnetization(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql
    AS $$
declare  
Begin
    if _params[3]=0 then
        return NULL;
    end if;
	RETURN ARRAY[_params[2]/_params[3]-_params[1]];
END

$$;


ALTER FUNCTION public.calc_magnetization(_params double precision[]) OWNER TO mena;

--
-- Name: calc_martensite_asc(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_martensite_asc(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$                                                                     
 import sympy as sy                                                                
 import math                                                                       
                                                                                   
 tem = _params[0]
 Mf = _params[1]
 Ms = _params[2]

 if  tem >= Ms:
        return [0]
 elif  tem <= Mf:
        return [1]
 else:
        return [1/2 * math.cos(math.pi * (tem-Mf)/(Ms-Mf))+1/2]                                                                                                       
 $$;


ALTER FUNCTION public.calc_martensite_asc(_params double precision[]) OWNER TO mena;

--
-- Name: calc_martensite_desc(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_martensite_desc(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$                                                                     
import math

tem = _params[0]
As = _params[1]
Af = _params[2]

if tem >= Af:
         return [0]
elif tem <= As:
         return [1]
else:
         return [1/2 * math.cos(math.pi * (tem-As)/(Af-As))+1/2]
$$;


ALTER FUNCTION public.calc_martensite_desc(_params double precision[]) OWNER TO mena;

--
-- Name: calc_max_block_pc(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_max_block_pc(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql STABLE
    AS $$                                              
 declare                                          
 	 d double precision = _params[1];
 	 sE double precision = _params[2];
 	 E_BFS double precision = _params[3];

 Begin 
       RETURN ARRAY[d/sE * E_BFS];
 
 END                                             
 $$;


ALTER FUNCTION public.calc_max_block_pc(_params double precision[]) OWNER TO mena;

--
-- Name: calc_max_block_stress_lh(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_max_block_stress_lh(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql
    AS $$                                                                         
 DECLARE                                                                               
         max_magn_stress double precision = _params[1];                                
         twin_stress double precision = _params[2];                                    
 BEGIN                                                            
         return ARRAY[max_magn_stress - twin_stress, max_magn_stress + twin_stress];   
 END                                                                                   
 $$;


ALTER FUNCTION public.calc_max_block_stress_lh(_params double precision[]) OWNER TO mena;

--
-- Name: calc_max_strain_mag(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_max_strain_mag(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql
    AS $$

 begin                                                                                                                      
     if _params[2]=0 then
         return NULL;
     end if;
     if 1.0-_params[1]/_params[2]<0 then
         return NULL;
     end if;
     return ARRAY[1.0-_params[1]/_params[2]];
 end;                                                                                                                       
 $$;


ALTER FUNCTION public.calc_max_strain_mag(_params double precision[]) OWNER TO mena;

--
-- Name: calc_neo_hook_minr(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_neo_hook_minr(ar double precision[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$
    
import math
import numpy as np
from scipy import optimize
from scipy import io
# [meas_time, mechanic_stress, strain]


_ar = np.array(ar)

lambda1 = np.array(_ar[:,2]+1)
sigma1 = np.array(_ar[:,1])


# Neo-Hooke
def NeoHooke_function(lam,Y):
    # Berechnet die technische Spannung anhand des Neo-Hooke-Modells
    sig = 1/lam*Y/3*(lam**2-1/(lam))
    return sig


# Definition des normierten Fehlers
def RNorm_function(X1,X2):
    R =  math.sqrt((X1 - X2)**2)/max(X1,0.01)
    return R

## Neo-Hooke
def NeoHooke_uni_fit_function_minR(sig,lam,Y):
    # gibt den mittleren normierten Fehler aus und wird zur Parameteridentifikation für den Elastizitätsmodul nach Neo-Hooke verwendet
    
    sig_model = NeoHooke_function(lam,Y)
    dist=0
    for i in range(0, len(lam)-1):
        dist =  RNorm_function(sig[i],sig_model[i])+dist
    
    dist=dist/len(lam)
    
    return dist

# Neo-Hooke-Modell
Y_NH_Start = 1e6 # Startwert
Y_NH_uni = optimize.fmin(lambda Y_NH:NeoHooke_uni_fit_function_minR(sigma1,lambda1,Y_NH), Y_NH_Start) # Y_NH_uni ist mein Ergbeniss des Fits
R_YNH=NeoHooke_uni_fit_function_minR(sigma1,lambda1,Y_NH_uni) # R ist der normierte Fehler als Zusatzinformtaion (Metadaten) für eine Aussage der Datenqualität 

return [Y_NH_uni[0]]
$$;


ALTER FUNCTION public.calc_neo_hook_minr(ar double precision[]) OWNER TO mena;

--
-- Name: calc_neo_hookean_initial_slope(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_neo_hookean_initial_slope(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql STABLE
    AS $$                                                                               
  declare                                                                                    
      sigma double precision = _params[1];                                                   
      eps double precision = _params[2];                                                     
                                                                                             
  Begin                                                                                      
          if ((1+eps)^2-1/(1+eps)) != 0 then                                                         
              return ARRAY[3*sigma/((1+eps)^2-1/(1+eps))];                                           
          end if;                                                                            
  END                                                                                        
  $$;


ALTER FUNCTION public.calc_neo_hookean_initial_slope(_params double precision[]) OWNER TO mena;

--
-- Name: calc_rel_perm(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_rel_perm(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql
    AS $$

declare  
Begin
  IF (_params[3]=0) OR (_params[2]=0) THEN
      return NULL;
  END IF;
  RETURN ARRAY[1+_params[1]/(_params[2]*_params[3])];
END

$$;


ALTER FUNCTION public.calc_rel_perm(_params double precision[]) OWNER TO mena;

--
-- Name: calc_yeoh_minr(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_yeoh_minr(ar double precision[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$
    
import math
import numpy as np
from scipy import optimize
from scipy import io
# [meas_time, mechanic_stress, strain]


_ar = np.array(ar)

lambda1 = np.array(_ar[:,2]+1)
sigma1 = np.array(_ar[:,1])



#Yeoh-Modell
def YeohUni_function(lam,C):
    # Berechnet die technische Spannung anhand des Yeoh-Modells
    I=(2/lam)+lam**2
    sig = 1/lam*(2*(lam**2)-2/(lam))*(C[0]+2*C[1]*(I-3)+3*C[2]*(I-3)**2)
    
    return sig


# Definition des normierten Fehlers
def RNorm_function(X1,X2):
    R =  math.sqrt((X1 - X2)**2)/max(X1,0.01)
    return R


# Yeoh-Modell
def Yeoh_uni_fit_function_minR(sig,lam,C):
    # gibt den mittleren normierten Fehler aus und wird zur Parameteridentifikation für die Yeoh-Parameter verwendet
    
    sig_model = YeohUni_function(lam,C)
    dist=0
    for i in range(0, len(lam)-1):
        dist =  RNorm_function(sig[i],sig_model[i])+dist
    
    dist=dist/len(lam)
    
    return dist

# Yeoh-Modell
C_Yeoh_Start = [0.3e6, 1e3, 1e3] # Startwert
C_Yeoh_uni = optimize.fmin(lambda C_Yeoh:Yeoh_uni_fit_function_minR(sigma1,lambda1,C_Yeoh), C_Yeoh_Start) # C_Yeoh_uni ist mein Ergbeniss des Fits
R_Yeoh=Yeoh_uni_fit_function_minR(sigma1,lambda1,C_Yeoh_uni) # R ist der normierte Fehler als Zusatzinformtaion (Metadaten) für eine Aussage der Datenqualität 

return [C_Yeoh_uni[0],C_Yeoh_uni[1],C_Yeoh_uni[2]]
$$;


ALTER FUNCTION public.calc_yeoh_minr(ar double precision[]) OWNER TO mena;

--
-- Name: calc_young_modulus_minr(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.calc_young_modulus_minr(ar double precision[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$
    
import math
import numpy as np
from scipy import optimize
from scipy import io
# [meas_time, mechanic_stress, strain]


_ar = np.array(ar)

lambda1 = np.array(_ar[:,2]+1)
sigma1 = np.array(_ar[:,1])


# linear Hooke
def Hooke_function(eps,Y):
    # Berechnet die technische Spannung anhand des linearen Hookschen Modells
    sig = Y*eps
    return sig

# Definition des normierten Fehlers
def RNorm_function(X1,X2):
    R =  math.sqrt((X1 - X2)**2)/max(X1,0.01)
    return R

# linear Hooke
def NeoHooke_fit_function_minR(sig,lam,Y):
    # gibt den mittleren normierten Fehler aus und wird zur Parameteridentifikation für den Elastizitätsmodul verwendet
    
    sig_model = Hooke_function(lam-1,Y)
    dist=0
    for i in range(0, len(lam)-1):
        dist =  RNorm_function(sig[i],sig_model[i])+dist
    
    dist=dist/len(lam)
    
    return dist


Y_Start = 1e6 # Startwert
Y = optimize.fmin(lambda Y:NeoHooke_fit_function_minR(sigma1,lambda1,Y), Y_Start) # Y ist mein Ergbeniss des Fits
R_Y=NeoHooke_fit_function_minR(sigma1,lambda1,Y) # R ist der normierte Fehler als Zusatzinformtaion (Metadaten) für eine Aussage der Datenqualität 

return [Y[0]]
$$;


ALTER FUNCTION public.calc_young_modulus_minr(ar double precision[]) OWNER TO mena;

--
-- Name: concr_param_function(integer, integer, integer); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.concr_param_function(start_step integer, end_step integer, _callversion integer) RETURNS TABLE(concr_param_id character varying, mod_meas_id character varying, param_id character varying, mat_sample_id character varying, meas_time integer, value double precision, mod_id character varying, step integer)
    LANGUAGE plpgsql PARALLEL RESTRICTED
    AS $$ DECLARE 
param_output param_values_type_output;

finished boolean;

_step integer;


initial_param record;

BEGIN

DELETE FROM temp_view;
DELETE FROM temp_view_results;
DELETE FROM temp_step_view;


FOR initial_param IN
SELECT
    concr_param.concr_param_id,
    concr_param.param_id,
    concr_param.mat_sample_id,
    concr_param.meas_time,
    concr_param.value
FROM
    concr_param
UNION
SELECT
    NULL :: character varying AS concr_param_id,
    constants.param_id,
    mat_sample.mat_sample_id,
    0 AS meas_time,
    constants.value
FROM
    constants,
    (
        SELECT
            sample.sample_id AS mat_sample_id
        FROM
            sample
        UNION
        SELECT
            material.mat_id AS mat_sample_id
        FROM
            material
    ) mat_sample
GROUP BY
    constants.param_id,
    mat_sample.mat_sample_id,
    constants.value
LOOP
    INSERT INTO
        temp_view(
            concr_param_id,
            mod_meas_id,
            param_id,
            mat_sample_id,
            meas_time,
            value,
            mod_id,
            step
        )
    SELECT
        initial_param.concr_param_id,
        NULL :: varchar AS mod_meas_id,
        initial_param.param_id,
        initial_param.mat_sample_id,
        initial_param.meas_time,
        initial_param.value,
        NULL :: varchar AS mod_id,
        0 AS step;

END LOOP;

_step = start_step + _callversion - 1;

finished = false;

WHILE (
    NOT finished
    AND _step < end_step
)

LOOP

finished = true;

_step = _step + 1;

DELETE FROM temp_view
    WHERE temp_view.concr_param_id IN (
        SELECT
            tsv.concr_param_id
        FROM 
            (SELECT
                tsv_inner.concr_param_id,
                tsv_inner.mod_meas_id,
                tsv_inner.param_id,
                tsv_inner.mat_sample_id,
                tsv_inner.meas_time,
                tsv_inner.value,
                tsv_inner.mod_id,
                k_1.matrix_indices
            FROM temp_view tsv_inner
            LEFT JOIN
                (
                    SELECT
                        DISTINCT pek1.concr_param_id,
                        ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                    FROM
                        concr_param_matrix pek1,
                        concr_param_matrix pek2
                    WHERE
                        pek1.concr_param_id = pek2.concr_param_id
                        AND pek1.index_id < pek2.index_id
                    UNION
                    SELECT
                        DISTINCT pek1.concr_param_id,
                        ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                    FROM
                        concr_param_matrix_additional pek1,
                        concr_param_matrix_additional pek2
                    WHERE
                        pek1.concr_param_id = pek2.concr_param_id
                        AND pek1.index_id < pek2.index_id
                ) k_1
                    ON tsv_inner.concr_param_id = k_1.concr_param_id
                WHERE tsv_inner.step < 1
            ) tsv
        LEFT JOIN 
            (SELECT
                tvr_inner.concr_param_id,
                tvr_inner.mod_meas_id,
                tvr_inner.param_id,
                tvr_inner.mat_sample_id,
                tvr_inner.meas_time,
                tvr_inner.value,
                tvr_inner.mod_id,
                k_1.matrix_indices
            FROM temp_view tvr_inner
            LEFT JOIN
                (
                    SELECT
                        DISTINCT pek1.concr_param_id,
                        ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                    FROM
                        concr_param_matrix pek1,
                        concr_param_matrix pek2
                    WHERE
                        pek1.concr_param_id = pek2.concr_param_id
                        AND pek1.index_id < pek2.index_id
                ) k_1
                    ON tvr_inner.concr_param_id = k_1.concr_param_id
                WHERE tvr_inner.step >= 1
            ) tv
                ON tsv.param_id = tv.param_id
                    AND tsv.mat_sample_id = tv.mat_sample_id
                    AND tsv.meas_time = tv.meas_time
                    AND (tsv.matrix_indices = tv.matrix_indices OR tsv.matrix_indices IS NULL)
            WHERE tv.param_id IS NOT NULL
    );

FOR param_output IN
SELECT
    (creating_result_set(array_agg(ROW(y.param_id, y.ar, y.mat_sample_id, y.mod_id, y.meas_time, y.mod_type)::param_values_type_input), _callversion)).concr_param_id AS concr_param_id,
    (creating_result_set(array_agg(ROW(y.param_id, y.ar, y.mat_sample_id, y.mod_id, y.meas_time, y.mod_type)::param_values_type_input), _callversion)).mod_meas_id AS mod_meas_id,
    (creating_result_set(array_agg(ROW(y.param_id, y.ar, y.mat_sample_id, y.mod_id, y.meas_time, y.mod_type)::param_values_type_input), _callversion)).param_id AS param_id,
    (creating_result_set(array_agg(ROW(y.param_id, y.ar, y.mat_sample_id, y.mod_id, y.meas_time, y.mod_type)::param_values_type_input), _callversion)).mat_sample_id AS mat_sample_id,
    (creating_result_set(array_agg(ROW(y.param_id, y.ar, y.mat_sample_id, y.mod_id, y.meas_time, y.mod_type)::param_values_type_input), _callversion)).meas_time AS meas_time,
    (creating_result_set(array_agg(ROW(y.param_id, y.ar, y.mat_sample_id, y.mod_id, y.meas_time, y.mod_type)::param_values_type_input), _callversion)).value AS value,
    (creating_result_set(array_agg(ROW(y.param_id, y.ar, y.mat_sample_id, y.mod_id, y.meas_time, y.mod_type)::param_values_type_input), _callversion)).mod_id AS mod_id
 
FROM
    (
        WITH RECURSIVE components_head_subs AS (
            SELECT
                cp1.concr_component_id AS parent,
                cp1.has_part AS child
            FROM
                component_parts cp1
                LEFT JOIN
                component_parts cp2
                	ON cp1.concr_component_id = cp2.has_part
                WHERE cp2.concr_component_id IS NULL
            UNION
            SELECT
                c.parent,
                cp.has_part AS child
            FROM
                components_head_subs c
                LEFT JOIN component_parts cp ON c.child = cp.concr_component_id
        ),
        mat_rec AS (
            
            SELECT
                m_param.param_id,
                m_param.value,
                m_param.mat_sample_id,
                ARRAY [m_param.mat_sample_id] AS mat_part_of,
                model.mod_id,
                model.mod_type,
                m_param.meas_time,
                m_input.max_in_number AS max_in_number,
                m_input.in_number,
                m_param.param_matrix_description,
                m_param.stepnr
            FROM
                (
                    SELECT
                        DISTINCT m_input_last.mod_id,
                        COALESCE(m_param_last.parent, last_round.mat_sample_id) AS mat_sample_id
                    FROM
                        (
                            SELECT
                                tlrv.param_id,
                                tlrv.mat_sample_id
                            FROM
                                temp_view tlrv
                            WHERE
                                tlrv.step = _step
                        ) last_round
                        LEFT JOIN (
                                SELECT
                                    DISTINCT
                                    param_piezo.param_id,
                                    param_piezo.param_matrix
                                FROM
                                    param_piezo
                             ) pp
                            ON pp.param_id = last_round.param_id
                        LEFT JOIN (
                            SELECT
                                modelinput.mod_id,
                                modelinput.param_id
                            FROM
                                modelinput
                            UNION
                            SELECT
                                model_condition_1.mod_id,
                                model_condition_1.param_id
                            FROM
                                model_condition model_condition_1
                            WHERE
                                model_condition_1.in_number < 0
                            UNION
                            SELECT
                                model_condition_2.mod_id,
                                model_condition_2.param_id
                            FROM
                                model_condition_compare model_condition_2
                            WHERE
                                model_condition_2.in_number < 0
                        ) m_input_last -- 					ON COALESCE(pp.param_id, last_round.param_id) = m_input_last.param_id
                        ON COALESCE(pp.param_matrix, last_round.param_id) = m_input_last.param_id
                        LEFT JOIN 
                            ( SELECT com.child, com.parent
                                FROM components_head_subs com
									WHERE com.child IS NOT NULL
                            ) m_param_last
                        ON m_param_last.child = last_round.mat_sample_id
                ) m_last_round
                LEFT JOIN (
                    SELECT
                        m_input_1.mod_id,
                        m_input_1.param_id,
                        count(m_input_1.in_number) OVER (PARTITION BY m_input_1.mod_id) AS max_in_number,
                        m_input_1.in_number
                    FROM
                        (
                            SELECT
                                modelinput.mod_id,
                                modelinput.param_id,
                                modelinput.in_number
                            FROM
                                modelinput 
                            UNION
                            SELECT
                                model_condition_1.mod_id,
                                model_condition_1.param_id,
                                model_condition_1.in_number
                            FROM
                                model_condition model_condition_1
                            WHERE
                                model_condition_1.in_number < 0
                            UNION
                            SELECT
                                model_condition_2.mod_id,
                                model_condition_2.param_id,
                                model_condition_2.in_number
                            FROM
                                model_condition_compare model_condition_2
                            WHERE
                                model_condition_2.in_number < 0
                        ) m_input_1
                ) m_input ON m_last_round.mod_id = m_input.mod_id
                LEFT JOIN components_head_subs ON components_head_subs.parent = m_last_round.mat_sample_id
                LEFT JOIN (
                    SELECT
                        CASE
                            WHEN k_1.matrix_indices IS NULL THEN m_param1.param_id
                            ELSE pp.param_matrix
                        END AS param_id,
                        m_param1.mat_sample_id,
                        m_param1.meas_time,
                        m_param1.value,
                        CASE
                            WHEN k_1.matrix_indices IS NOT NULL THEN ARRAY [m_param1.param_id] || k_1.matrix_indices
                            ELSE NULL :: character varying []
                        END AS param_matrix_description,
                        m_param1.stepnr
                    FROM
                        (
                            SELECT
                                param_values_1.concr_param_id,
                                param_values_1.param_id,
                                param_values_1.mat_sample_id,
                                param_values_1.meas_time,
                                param_values_1.value,
                                param_values_1.step AS stepnr
                            FROM
                                temp_view param_values_1
                            WHERE
                                param_values_1.step <= end_step
                       ) m_param1
                        LEFT JOIN (
                            SELECT
                                DISTINCT
                                param_piezo.param_id,
                                param_piezo.param_matrix
                            FROM
                                param_piezo
                        ) pp ON pp.param_id = m_param1.param_id
                        LEFT JOIN (
                            SELECT
                                DISTINCT pek1.concr_param_id,
                                ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                            FROM
                                concr_param_matrix pek1,
                                concr_param_matrix pek2
                            WHERE
                                pek1.concr_param_id = pek2.concr_param_id
                                AND pek1.index_id < pek2.index_id
                            UNION
                            SELECT
                                DISTINCT pek1.concr_param_id,
                                ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                            FROM
                                concr_param_matrix_additional pek1,
                                concr_param_matrix_additional pek2
                            WHERE
                                pek1.concr_param_id = pek2.concr_param_id
                                AND pek1.index_id < pek2.index_id
                        ) k_1 ON m_param1.concr_param_id = k_1.concr_param_id
                ) m_param ON (
                    m_input.param_id :: text = m_param.param_id :: text
                    OR m_input.param_id :: text = m_param.param_matrix_description [1] :: text
                )
                AND COALESCE(components_head_subs.child, m_last_round.mat_sample_id) = m_param.mat_sample_id
                LEFT JOIN model ON model.mod_id = m_input.mod_id
                LEFT JOIN (
                    SELECT
                        DISTINCT mie1.param_id,
                        mie1.mod_id,
                        ARRAY [mie1.param_id::character varying, mie1.index_id::character varying,
                                mie1.value::character varying,
                                mie2.index_id::character varying,
                                mie2.value::character varying] AS input_indices
                    FROM
                        modelinput_matrix mie1,
                        modelinput_matrix mie2
                    WHERE
                        mie1.param_id :: text = mie2.param_id :: text
                        AND mie1.mod_id :: text = mie2.mod_id :: text
                        AND mie1.index_id < mie2.index_id
                ) k_2 ON m_param.param_matrix_description :: text = k_2.input_indices :: text
                AND m_input.mod_id = k_2.mod_id
            WHERE
                (
                    m_input.mod_id not IN (
                        Select
                            modelinput_matrix.mod_id
                        From
                            modelinput_matrix
                    )
                    or m_param.param_matrix_description :: text = k_2.input_indices :: text
                    or m_param.param_matrix_description :: text is NULL
                    
                )
                AND NOT m_param.value IS NULL
            UNION
      
            SELECT
             
                mat_rec_1.param_id,
                mat_rec_1.value,
                component_parts.concr_component_id AS mat_sample_id,
                ARRAY [component_parts.concr_component_id] || mat_rec_1.mat_part_of AS mat_part_of,
                mat_rec_1.mod_id,
                mat_rec_1.mod_type,
                mat_rec_1.meas_time,
                mat_rec_1.max_in_number,
                mat_rec_1.in_number,
                mat_rec_1.param_matrix_description,
                mat_rec_1.stepnr
            FROM
                mat_rec mat_rec_1,
                component_parts
            WHERE
                mat_rec_1.mat_sample_id = component_parts.has_part
                AND mat_rec_1.meas_time = 0
        ) 
        SELECT
           
            CASE
                WHEN NOT mat_rec_all.mod_type = 'matrixoperation' THEN model_out.param_id
                ELSE ARRAY [get_matrix_output_param_id(
                        k.matrix_indices,
                        modeloutput_pek_bits."bit",
                        model_out.param_id[1]
            ) ]
    END AS param_id,
    
    CASE
        WHEN NOT mat_rec_all.mod_type = 'matrixoperation' THEN array_agg(
            mat_rec_all.value
            ORDER BY
                mat_rec_all.in_number
        ) FILTER (
            WHERE
                mat_rec_all.in_number > 0
        ) :: character varying []
        ELSE array_agg(
            (
                (
                    ARRAY [mat_rec_all.param_id] || mat_rec_all.param_matrix_description [2:]
                ) || ARRAY [mat_rec_all.value::character varying]
            ) || k.matrix_indices
            ORDER BY
                mat_rec_all.in_number
        ) FILTER (
            WHERE
                mat_rec_all.in_number > 0
        )
    END AS ar,
    
     ARRAY [mat_rec_all.mat_sample_id] AS mat_sample_id,
    mat_rec_all.mod_id,
    mat_rec_all.mod_type,
    mat_rec_all.meas_time,
    array_agg(
        mat_rec_all.value
        ORDER BY
            mat_rec_all.in_number
    ) FILTER (
        WHERE
            mat_rec_all.in_number <= -100
    ) :: character varying [] AS groupbyfilter,
    max(mat_rec_all.stepnr) AS stepnr
FROM
    (
        SELECT
            mat_rec_all_1.param_id,
            mat_rec_all_1.value,
            mat_rec_all_1.mat_sample_id,
            mat_rec_all_1.mat_part_of,
            mat_rec_all_1.mod_id,
            mat_rec_all_1.mod_type,
            mat_rec_all_1.meas_time,
            mat_rec_all_1.max_in_number,
            mat_rec_all_1.in_number,
            mat_rec_all_1.param_matrix_description,
            mat_rec_all_1.stepnr
        FROM
            mat_rec mat_rec_all_1 
        UNION
        SELECT
            mat_rec_0.param_id,
            mat_rec_0.value,
            mat_rec_0.mat_sample_id,
            mat_rec_0.mat_part_of,
            mat_rec_0.mod_id,
            mat_rec_0.mod_type,
            mat_rec_1.meas_time,
            mat_rec_0.max_in_number,
            mat_rec_0.in_number,
            mat_rec_0.param_matrix_description,
            mat_rec_0.stepnr
        FROM
            mat_rec mat_rec_0,
            mat_rec mat_rec_1
        WHERE
            mat_rec_0.meas_time = 0
            AND mat_rec_0.mat_sample_id = mat_rec_1.mat_sample_id
            AND mat_rec_0.mod_id = mat_rec_1.mod_id
    ) mat_rec_all
    LEFT JOIN model_condition ON mat_rec_all.mod_id = model_condition.mod_id
    AND mat_rec_all.param_id = model_condition.param_id
     LEFT JOIN model_condition_compare ON mat_rec_all.mod_id = model_condition_compare.mod_id
    AND mat_rec_all.param_id = model_condition_compare.param_id AND  mat_rec_all.in_number = model_condition_compare.in_number
    LEFT JOIN (
        SELECT
            DISTINCT pek1.param_id,
            ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
        FROM
            param_piezo pek1,
            param_piezo pek2
        WHERE
            pek1.param_id = pek2.param_id
            AND pek1.index_id < pek2.index_id
    ) k ON k.param_id = mat_rec_all.param_matrix_description [1]-- and mat_rec_all.mod_type = 'matrixoperation'
    LEFT JOIN (
        SELECT
            array_agg(
                modeloutput.param_id
                ORDER BY
                    out_number
            ) AS param_id,
            modeloutput.mod_id
        FROM
            modeloutput
        GROUP BY
            modeloutput.mod_id
    ) model_out ON mat_rec_all.mod_id = model_out.mod_id
    LEFT JOIN modeloutput_pek_bits ON model_out.mod_id = modeloutput_pek_bits.mod_id 
 
WHERE
    model_condition.value IS NULL
    OR test_condition(
        mat_rec_all.value,
        model_condition.value :: double precision,
        model_condition.operator
    ) 
GROUP BY
    mat_rec_all.mod_type,
    ARRAY [mat_rec_all.mat_sample_id],
    model_out.param_id,
    mat_rec_all.mod_id,
    mat_rec_all.meas_time,
    k.matrix_indices,
    modeloutput_pek_bits."bit"

HAVING

        ((
        count(mat_rec_all.value) = max(mat_rec_all.max_in_number)
        AND is_distinct(
            array_agg(mat_rec_all.in_number) FILTER (where mat_rec_all.in_number > 0) :: double precision []
        )
        OR mat_rec_all.mod_type = 'linear_interpolation'
        AND (max(mat_rec_all.in_number) = 1 or count(mat_rec_all.value) = max(mat_rec_all.max_in_number) -1)
        OR mat_rec_all.mod_type = 'matrixoperation'
        AND count(mat_rec_all.value) > 0
        AND matrix_filter_agg(
            (
                mat_rec_all.stepnr,
                mat_rec_all.param_matrix_description
            ) :: matrix_filter_agg_type
        ) = _step
    )
    AND (size_1_or_2nd_diff(mat_rec_all.mat_part_of)
        OR ARRAY[(max(mat_rec_all.mat_part_of))[1], mat_rec_all.mod_id] in 
                (Select ARRAY[concr_component_id, material_condition.mod_id] 
                from concr_component, material_condition 
                where component_id = mat_type )
    ))
   AND 
    
    (CASE WHEN max(model_condition_compare.mod_id) is NULL THEN true 
        ELSE 
     test_compare_condition(array_agg( mat_rec_all.value ORDER BY model_condition_compare.in_number, model_condition_compare.compare_position ) FILTER (WHERE
        model_condition_compare.in_number < 0 ), array_agg(model_condition_compare.operator ORDER BY model_condition_compare.in_number, model_condition_compare.compare_position) FILTER (WHERE
        model_condition_compare.in_number < 0))
         END)

   
) y
LEFT JOIN (
    SELECT
        material_condition.mod_id,
        material.mat_id AS mat_sample_id
    FROM
        material_condition
        left join material on (
            material_condition.mat_type :: text = material.mat_type :: text
        )
    UNION
    SELECT
        material_condition.mod_id,
        sample_id AS mat_sample_id
    FROM
        material_condition
        left join (
            Select
                sample_id,
                mat_type
            from
                material,
                sample
            where
                sample.mat_id :: text = material.mat_id :: text
        ) m on (
            material_condition.mat_type :: text = m.mat_type :: text
        )
    UNION
    Select
        material_condition.mod_id,
        concr_component.concr_component_id AS mat_sample_id
    from
        material_condition
        left join concr_component on (
            material_condition.mat_type = concr_component.component_id
        )
) mat_condition ON mat_condition.mod_id :: text = y.mod_id :: text
WHERE
    (
        mat_condition.mod_id is null
        or mat_condition.mat_sample_id = y.mat_sample_id [1]
    ) 
GROUP BY
    y.mod_id,
    y.mod_type,
    y.param_id,
    (y.mat_sample_id) [1],
    y.groupbyfilter,
   
    CASE
        WHEN NOT y.mod_type LIKE 'linear_interpolation'
        AND NOT y.mod_type LIKE 'agg%'
        AND NOT y.mod_type LIKE 'complex_sum' THEN y.meas_time :: character varying
        ELSE NULL :: character varying
    END
HAVING
    max(y.stepnr) = _step
    AND size_1_or_2nd_diff(y.mat_sample_id) 
        
    LOOP finished = false;

INSERT INTO
    temp_step_view(
        concr_param_id,
        mod_meas_id,
        param_id,
        mat_sample_id,
        meas_time,
        value,
        mod_id
    )
SELECT
    param_output.concr_param_id,
    param_output.mod_meas_id,
    param_output.param_id,
    param_output.mat_sample_id,
    param_output.meas_time,
    param_output.value,
    param_output.mod_id;

END LOOP;


-- Write Data for the Next Round (eliminated duplicates exkl. mod_id)
INSERT INTO
    temp_view(
        concr_param_id,
        mod_meas_id,
        param_id,
        mat_sample_id,
        meas_time,
        value,
        mod_id,
        step
    )
SELECT
		(array_agg(tsv.concr_param_id ORDER BY tsv.ctid))[1] AS concr_param_id,
		(array_agg(tsv.mod_meas_id ORDER BY tsv.ctid))[1] AS mod_meas_id,
        tsv.param_id,
        tsv.mat_sample_id,
        tsv.meas_time,
		(array_agg(tsv.value ORDER BY tsv.ctid))[1] AS value,
		(array_agg(tsv.mod_id ORDER BY tsv.ctid))[1] AS mod_id,
        _step + 1 AS step
    FROM 
        (SELECT
            tsv_inner.concr_param_id,
            tsv_inner.mod_meas_id,
            tsv_inner.param_id,
            tsv_inner.mat_sample_id,
            tsv_inner.meas_time,
            tsv_inner.value,
            tsv_inner.mod_id,
            tsv_inner.ctid,
            k_1.matrix_indices
        FROM temp_step_view tsv_inner
        LEFT JOIN
            (
                SELECT
                    DISTINCT pek1.concr_param_id,
                    ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                FROM
                    concr_param_matrix pek1,
                    concr_param_matrix pek2
                WHERE
                    pek1.concr_param_id = pek2.concr_param_id
                    AND pek1.index_id < pek2.index_id
                UNION
                SELECT
                    DISTINCT pek1.concr_param_id,
                    ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                FROM
                    concr_param_matrix_additional pek1,
                    concr_param_matrix_additional pek2
                WHERE
                    pek1.concr_param_id = pek2.concr_param_id
                    AND pek1.index_id < pek2.index_id
            ) k_1
                ON tsv_inner.concr_param_id = k_1.concr_param_id
        ) tsv
    LEFT JOIN 
        (SELECT
            tvr_inner.concr_param_id,
            tvr_inner.mod_meas_id,
            tvr_inner.param_id,
            tvr_inner.mat_sample_id,
            tvr_inner.meas_time,
            tvr_inner.value,
            tvr_inner.mod_id,
            k_1.matrix_indices
        FROM temp_view_results tvr_inner
        LEFT JOIN
            (
                SELECT
                    DISTINCT pek1.concr_param_id,
                    ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                FROM
                    concr_param_matrix pek1,
                    concr_param_matrix pek2
                WHERE
                    pek1.concr_param_id = pek2.concr_param_id
                    AND pek1.index_id < pek2.index_id
                UNION
                SELECT
                    DISTINCT pek1.concr_param_id,
                    ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                FROM
                    concr_param_matrix_additional pek1,
                    concr_param_matrix_additional pek2
                WHERE
                    pek1.concr_param_id = pek2.concr_param_id
                    AND pek1.index_id < pek2.index_id
            ) k_1
                ON tvr_inner.concr_param_id = k_1.concr_param_id
        ) tv
            ON tsv.param_id = tv.param_id
                AND tsv.mat_sample_id = tv.mat_sample_id
                AND tsv.meas_time = tv.meas_time
                AND (tsv.matrix_indices = tv.matrix_indices OR tsv.matrix_indices IS NULL)
        WHERE tv.param_id IS NULL
        GROUP BY
			tsv.param_id,
			tsv.mat_sample_id,
			tsv.meas_time,
            tsv.matrix_indices;


INSERT INTO
    temp_view_results(
        concr_param_id,
        mod_meas_id,
        param_id,
        mat_sample_id,
        meas_time,
        value,
        mod_id,
        step
    )
SELECT
        tsv.concr_param_id,
        tsv.mod_meas_id,
        tsv.param_id,
        tsv.mat_sample_id,
        tsv.meas_time,
        tsv.value,
        tsv.mod_id,
        _step + 1 AS step
    FROM 
        (SELECT
            tsv_inner.concr_param_id,
            tsv_inner.mod_meas_id,
            tsv_inner.param_id,
            tsv_inner.mat_sample_id,
            tsv_inner.meas_time,
            tsv_inner.value,
            tsv_inner.mod_id,
            k_1.matrix_indices
        FROM temp_step_view tsv_inner
        LEFT JOIN
            (
                SELECT
                    DISTINCT pek1.concr_param_id,
                    ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                FROM
                    concr_param_matrix pek1,
                    concr_param_matrix pek2
                WHERE
                    pek1.concr_param_id = pek2.concr_param_id
                    AND pek1.index_id < pek2.index_id
                UNION
                SELECT
                    DISTINCT pek1.concr_param_id,
                    ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                FROM
                    concr_param_matrix_additional pek1,
                    concr_param_matrix_additional pek2
                WHERE
                    pek1.concr_param_id = pek2.concr_param_id
                    AND pek1.index_id < pek2.index_id
            ) k_1
                ON tsv_inner.concr_param_id = k_1.concr_param_id
        ) tsv
    LEFT JOIN 
        (SELECT
            tvr_inner.concr_param_id,
            tvr_inner.mod_meas_id,
            tvr_inner.param_id,
            tvr_inner.mat_sample_id,
            tvr_inner.meas_time,
            tvr_inner.value,
            tvr_inner.mod_id,
            k_1.matrix_indices
        FROM temp_view_results tvr_inner
        LEFT JOIN
            (
                SELECT
                    DISTINCT pek1.concr_param_id,
                    ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                FROM
                    concr_param_matrix pek1,
                    concr_param_matrix pek2
                WHERE
                    pek1.concr_param_id = pek2.concr_param_id
                    AND pek1.index_id < pek2.index_id
                UNION
                SELECT
                    DISTINCT pek1.concr_param_id,
                    ARRAY [pek1.index_id::character varying, pek1.value::character varying, pek2.index_id::character varying, pek2.value::character varying] AS matrix_indices
                FROM
                    concr_param_matrix_additional pek1,
                    concr_param_matrix_additional pek2
                WHERE
                    pek1.concr_param_id = pek2.concr_param_id
                    AND pek1.index_id < pek2.index_id
            ) k_1
                ON tvr_inner.concr_param_id = k_1.concr_param_id
        ) tv
            ON tsv.param_id = tv.param_id
                AND tsv.mat_sample_id = tv.mat_sample_id
                AND tsv.meas_time = tv.meas_time
                AND tsv.mod_id = tv.mod_id
                AND (tsv.matrix_indices = tv.matrix_indices OR tsv.matrix_indices IS NULL)
        WHERE tv.param_id IS NULL;

DELETE FROM temp_step_view;



END LOOP;

RETURN QUERY
SELECT
    tv.concr_param_id,
    tv.mod_meas_id,
    tv.param_id,
    tv.mat_sample_id,
    tv.meas_time,
    tv.value,
    tv.mod_id,
    tv.step
FROM
    temp_view_results tv;
END $$;


ALTER FUNCTION public.concr_param_function(start_step integer, end_step integer, _callversion integer) OWNER TO mena;

--
-- Name: creating_result_set(public.param_values_type_input[], integer); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.creating_result_set(_params public.param_values_type_input[], _callversion integer) RETURNS SETOF public.param_values_type_output
    LANGUAGE plpgsql
    AS $$
DECLARE
	row_content param_values_type_output;
	row_content_matrix param_values_type_output_matrix;
	row_content_indices param_values_type_output_matrix;
	row_content_meas_time param_meas_time_type_creating_intern;
	_mod_type varchar;
	_mod_id varchar;
	_param_id varchar;
	_card int;
	row_content_aggregat aggregat_type_creating_intern;

BEGIN
	SELECT mod_type from unnest(_params) INTO _mod_type;
	SELECT mod_id from unnest(_params) INTO _mod_id;

	IF _mod_type = 'simple_calculation' THEN
		FOR row_content_meas_time in
			SELECT 'cp_mod_'||param_id[1] ||'_'||model.mod_id||'_'|| mat_sample_id[1] ||'_'||meas_time as concr_param_id_pre,
			'mod_'||model.mod_id  ||'_'||param_id[1] ||'_'|| mat_sample_id[1]||'_'||meas_time  as mod_meas_id,
			param_id,
			mat_sample_id[1] AS mat_sample_id,
			meas_time,
			simple_function(function_name, cast(ar as double precision[])) as value,
			model.mod_id
		FROM	unnest(_params) p,
			model
		WHERE 	model.mod_id = p.mod_id

		LOOP
			
			If row_content_meas_time.value IS NOT NULL THEN

				FOR i IN 1 .. array_length((row_content_meas_time.param_ids),1)
				
				LOOP
				
					If row_content_meas_time.value[i] IS NOT NULL THEN
						row_content.concr_param_id = row_content_meas_time.concr_param_id_pre || '_' || (row_content_meas_time.param_ids)[i];
						row_content.param_id = (row_content_meas_time.param_ids)[i];
						row_content.mod_meas_id = row_content_meas_time.mod_meas_id_pre  || '_' || (row_content_meas_time.param_ids)[i];
						row_content.mat_sample_id = row_content_meas_time.mat_sample_id;
						row_content.meas_time = row_content_meas_time.meas_time;
						row_content.value = row_content_meas_time.value[i];
						row_content.mod_id = row_content_meas_time.mod_id;
						RETURN NEXT row_content;
					END IF;

				END LOOP;
			ELSE
				RETURN;
			END IF;
		END LOOP;
	ELSIF _mod_type = 'complex_sum' THEN
		
		FOR row_content_matrix in
                SELECT 'cp_'||  'mod_'||'_'||param_id[1] ||'_'||p.mod_id||'_'|| mat_sample_id[1] ||'_'||max(ar[1])  as concr_param_id,
                         'mod_'||p.mod_id  ||'_'||param_id[1] ||'_'|| mat_sample_id[1] ||'_'||max(ar[1])  as mod_meas_id,
                        param_id[1] AS param_id,
                        mat_sample_id[1] AS mat_sample_id,
                        sum(meas_time)+1000  as meas_time,
                      CASE
			WHEN count(ar)=2 THEN
				 integral_function(function_name, array_agg(cast(ar as double precision[])))
        		ELSE
				NULL
			END value,
	                p.mod_id
                FROM 	unnest(_params) p, model
		WHERE p.mod_id = model.mod_id
		GROUP BY param_id, p.mod_id, mat_sample_id[1], function_name, ar[1]
	        LOOP
        	       If row_content_matrix.value != '{}' THEN
				
                                SELECT modelinput.param_id FROM modelinput, unnest(_params) p  where in_number = 1 and modelinput.mod_id =p.mod_id INTO _param_id;
				row_content.param_id = _param_id;
				row_content.concr_param_id = row_content_matrix.concr_param_id || _param_id;
                                row_content.mat_sample_id = row_content_matrix.mat_sample_id;
                                row_content.meas_time = row_content_matrix.meas_time;
                                row_content.value = row_content_matrix.value[1];
	                        row_content.mod_id = row_content_matrix.mod_id;
                       		RETURN NEXT row_content;

                                row_content.concr_param_id = row_content_matrix.concr_param_id;
                                row_content.param_id = row_content_matrix.param_id;
                                row_content.mat_sample_id = row_content_matrix.mat_sample_id;
                                row_content.meas_time = row_content_matrix.meas_time;
                                row_content.value = row_content_matrix.value[2];
                                row_content.mod_id = row_content_matrix.mod_id;

                                RETURN NEXT row_content;

           
                	END IF;
        	END LOOP;
	ELSIF _mod_type = 'agg_integral' THEN

	FOR row_content_aggregat in
		SELECT
			
			 'cp_'||  'mod_'|| mod_id ||'_' || mat_sample_id[1] || '_' as concr_param_id_pre,
			 'mod_'|| mod_id  ||'_' || mat_sample_id[1] ||'_' as mod_meas_id_pre,
			
			param_id AS param_ids,
			mat_sample_id[1] AS mat_sample_id,
			
			value,
			mod_id
		FROM (	SELECT unnest_nd_1d(integral_function(max(function_name), array_agg(meas_time || cast(ar as double precision[]) order by meas_time))) as value
			FROM unnest(_params) x, model WHERE model.mod_id = x.mod_id ) calc, unnest(_params) params
		WHERE calc.value[1] = params.meas_time
		LOOP
				
				IF row_content_aggregat.value != '{}'  THEN
					
						FOR i IN 1 .. array_length((row_content_aggregat.param_ids),1) LOOP
								IF row_content_aggregat.value[i+2] IS NOT NULL THEN
										row_content.concr_param_id = row_content_aggregat.concr_param_id_pre || (row_content_aggregat.param_ids)[i] || '_' || row_content_aggregat.value[2];
										row_content.param_id = (row_content_aggregat.param_ids)[i];
										row_content.mod_meas_id = row_content_aggregat.mod_meas_id_pre  || (row_content_aggregat.param_ids)[i];
										row_content.mat_sample_id = row_content_aggregat.mat_sample_id;
										row_content.meas_time = row_content_aggregat.value[2];
										row_content.value = row_content_aggregat.value[i+2];
										row_content.mod_id = row_content_aggregat.mod_id;
										RETURN NEXT row_content;
								
							END IF;
				  END LOOP;
				
				END IF;
		END LOOP;

	ELSIF _mod_type = 'aggregate' THEN

	FOR row_content_aggregat in
		SELECT
			
			 'cp_'||  'mod_'|| mod_id ||'_' || mat_sample_id[1]  as concr_param_id_pre,
			 'mod_'|| mod_id  ||'_' || mat_sample_id[1]  as mod_meas_id_pre,
			
			param_id AS param_ids,
			mat_sample_id[1] AS mat_sample_id,
			
			value,
			mod_id
		FROM (	SELECT simple_function(max(function_name), array_agg(meas_time || cast(ar as double precision[]) order by meas_time)) as value, max(param_id) as param_id, max(mat_sample_id) as mat_sample_id, max(model.mod_id) as mod_id
			FROM unnest(_params) x, model WHERE model.mod_id = x.mod_id ) calc
		LOOP
				
				IF row_content_aggregat.value != '{}'  THEN

								FOR i IN 1 .. array_length((row_content_aggregat.param_ids),1) LOOP
										row_content.concr_param_id = row_content_aggregat.concr_param_id_pre || (row_content_aggregat.param_ids)[i];
										row_content.param_id = (row_content_aggregat.param_ids)[i];
										row_content.mod_meas_id = row_content_aggregat.mod_meas_id_pre  || (row_content_aggregat.param_ids)[i];
										row_content.mat_sample_id = row_content_aggregat.mat_sample_id;
										row_content.meas_time = 0;
										row_content.value = row_content_aggregat.value[i];
										row_content.mod_id = row_content_aggregat.mod_id;
										RETURN NEXT row_content;
							END LOOP;
				
				END IF;
		END LOOP;
	
	ELSIF _mod_type = 'matrixoperation' THEN
		FOR row_content_matrix in
			SELECT 'cp_'||  'mod_'||param_id[1] ||'_'||model.mod_id||'_'|| mat_sample_id[1] ||'_'||meas_time as concr_param_id,
			'mod_'||model.mod_id  ||'_'||param_id[1] ||'_'|| mat_sample_id[1] ||'_'||meas_time  as mod_meas_id,
			param_id[1] AS param_id,
			mat_sample_id[1] AS mat_sample_id,
			meas_time,
			matrix_function(function_name, ar) as value,
			model.mod_id
		FROM
			unnest(_params) p,
			model
		WHERE 	model.mod_id = p.mod_id

		LOOP
			If row_content_matrix.value != '{}' THEN
				
				FOR i IN 1 .. array_length(row_content_matrix.value, 1) LOOP
					row_content_indices = row_content_matrix;
					row_content_indices.concr_param_id = row_content_matrix.concr_param_id ||'_'|| i;
					row_content_indices.value = ARRAY[row_content_matrix.value[i][5]];
				
					row_content = NULL;
					row_content.concr_param_id = row_content_indices.concr_param_id;
					row_content.param_id = row_content_indices.param_id;
					row_content.mat_sample_id = row_content_indices.mat_sample_id;
					row_content.meas_time = row_content_indices.meas_time;
					row_content.value = row_content_indices.value[1];
					row_content.mod_id = row_content_indices.mod_id;
					row_content.mod_meas_id = row_content_indices.mod_meas_id;

					INSERT INTO concr_param_matrix_additional VALUES (row_content_indices.concr_param_id, row_content_matrix.value[i][1], row_content_matrix.value[i][2]) ON CONFLICT ON CONSTRAINT concr_param_matrix_additional_pkey DO UPDATE set value = row_content_matrix.value[i][2];
					INSERT INTO concr_param_matrix_additional VALUES (row_content_indices.concr_param_id, row_content_matrix.value[i][3], row_content_matrix.value[i][4]) ON CONFLICT ON CONSTRAINT concr_param_matrix_additional_pkey DO UPDATE set value = row_content_matrix.value[i][4];
					RETURN NEXT row_content;
				END LOOP;
			ELSE
				RETURN;
			END IF;
		END LOOP;

		
    ELSIF _mod_type = 'aggregate_simple' THEN
		FOR row_content in
			SELECT 'cp_'||  'mod_'||model.mod_id  ||'_'||param_id[1] ||'_'|| mat_sample_id[1] ||'_'||0 as concr_param_id,
			 'mod_'||model.mod_id  ||'_'||param_id[1] ||'_'|| mat_sample_id[1] as mod_meas_id,
			param_id[1] AS param_id,
			mat_sample_id[1] AS mat_sample_id,
			0 as meas_time,
			aggregat_simple(function_name, array_agg(cast(ar[1] as double precision))) as value,
			model.mod_id
		FROM	unnest(_params) cp,
			model
		where model.mod_id = cp.mod_id
		
		group by mat_sample_id[1], param_id, model.mod_id, function_name
		LOOP
			If row_content.value IS NOT NULL THEN
				RETURN NEXT row_content;
			END IF;
		END LOOP;
    
    ELSIF _mod_type = 'linear_interpolation' THEN
        SELECT count(*) from unnest(_params) where cardinality(ar) = 3 into _card;
		IF _card < 2 then
			return;
        ELSE
		FOR row_content in
			SELECT 'cp_mod_'||function_name||'_'||param_id[1] ||'_'||model.mod_id||'_'|| mat_sample_id[1] ||'_0'  as concr_param_id,
			'mod_'||model.mod_id  ||'_'||param_id[1] ||'_'|| mat_sample_id[1] ||'_0'  as mod_meas_id,
			param_id[1] AS param_id,
			mat_sample_id[1] AS mat_sample_id,
			0 as meas_time,
			interpol_function(function_name, array_agg(cast((cast(ar[1] as double precision), cast(ar[2] as double precision)) as interpol_type)) FILTER (where ar[1] is not null) ) as value,
			model.mod_id
			FROM
				(Select param_id,
	                    ar[2:3] as ar,
                        mat_sample_id,
                        mod_id
                FROM
                unnest(_params) 
                UNION ALL
                Select distinct  param_id,
	                    ARRAY[ar[1]]  as ar,
                        mat_sample_id,
                        mod_id
                        from 
                        unnest(_params) ) p,
			    model
		    WHERE 	model.mod_id = p.mod_id
            group by function_name, param_id, model.mod_id, mat_sample_id
		LOOP
			If row_content.value IS NOT NULL THEN
				RETURN NEXT row_content;
			END IF;
		END LOOP;
		END IF;
	END IF;
	RETURN;
	END
$$;


ALTER FUNCTION public.creating_result_set(_params public.param_values_type_input[], _callversion integer) OWNER TO mena;

--
-- Name: function_concr_curve_params(); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.function_concr_curve_params() RETURNS TABLE(id_tuple text, id_curve text, curve_id character varying, mat_sample_id character varying, meas_time integer)
    LANGUAGE sql IMMUTABLE PARALLEL RESTRICTED
    AS $$
SELECT curve_params.id_curve||'_'||curve_params.curve_id||'_'||curve_params.mat_sample_id||'_'||unnest(curve_params.meas_time_array) AS id_tuple,
   curve_params.id_curve,
   curve_params.curve_id,
   curve_params.mat_sample_id,
   unnest(curve_params.meas_time_array) AS meas_time
   FROM ( SELECT curve_n.curve_id||'_'||curve_n.mat_sample_id||'_'||array_to_string(array_agg(curve_n.meas_time),'_') AS id_curve,
            curve_n.curve_id,
            curve_n.mat_sample_id,
            array_agg(curve_n.meas_time) AS meas_time_array
         FROM ( SELECT p.mat_sample_id,
                  p.meas_time,
                  array_agg(p.param_id ORDER BY p.param_id) AS params_given,
                  array_agg(p.value ORDER BY p.param_id) AS params_value,
                  curve.curve_id
                  FROM concr_parameter p,
                  characteristic_curve curve
                  GROUP BY curve.curve_id, p.mat_sample_id, p.meas_time) curve_n
            JOIN ( SELECT curve_params_2.curve_id,
                  array_agg(curve_params_2.param_id) AS params_needed
                  FROM curve_params curve_params_2
                  GROUP BY curve_params_2.curve_id) curve_params_1 ON curve_n.curve_id::text = curve_params_1.curve_id::text
            LEFT JOIN curve_condition ON curve_condition.curve_id::text = curve_params_1.curve_id::text
         WHERE array_contained(curve_params_1.params_needed, curve_n.params_given, 2) AND (curve_condition.curve_id IS NULL OR (curve_condition.param_id::text = ANY (curve_n.params_given::text[])) AND (curve_condition.in_number <= '-100'::integer OR curve_n.params_value[array_position(curve_n.params_given, curve_condition.param_id)] = curve_condition.value::double precision))
         GROUP BY curve_n.curve_id, curve_n.mat_sample_id, (
               CASE
                  WHEN curve_condition.param_id IS NOT NULL THEN curve_n.params_value[array_position(curve_n.params_given, curve_condition.param_id)]
                  ELSE NULL::double precision
               END)
         HAVING count(*) > 1) curve_params;
$$;


ALTER FUNCTION public.function_concr_curve_params() OWNER TO mena;

--
-- Name: function_concr_curve_params_temp(); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.function_concr_curve_params_temp() RETURNS TABLE(concr_tuple_id text, id_tuple text, concr_curve_id text, concr_param_id1 text, concr_param_id2 text)
    LANGUAGE sql IMMUTABLE PARALLEL RESTRICTED
    AS $$
 SELECT DISTINCT ON ((concr_curve_params.id_tuple || concr_curve.concr_curve_id)) concr_curve_params.id_tuple || concr_curve.concr_curve_id AS concr_tuple_id,
    concr_curve_params.id_tuple,
    concr_curve.concr_curve_id,
    max(concr_parameter.concr_param_id::text) FILTER (WHERE concr_curve.in_number = 1) AS concr_param_id1,
    max(concr_parameter.concr_param_id::text) FILTER (WHERE concr_curve.in_number = 2) AS concr_param_id2
   FROM concr_curve,
    concr_curve_params,
    concr_parameter
  WHERE concr_curve.id_curve = concr_curve_params.id_curve AND concr_parameter.mat_sample_id::text = concr_curve_params.mat_sample_id::text AND concr_parameter.meas_time = concr_curve_params.meas_time AND concr_parameter.param_id::text = concr_curve.param_id::text
  GROUP BY concr_curve_params.id_tuple, concr_curve.concr_curve_id;
$$;


ALTER FUNCTION public.function_concr_curve_params_temp() OWNER TO mena;

--
-- Name: get_matrix_output_param_id(character varying[], character varying, character varying); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.get_matrix_output_param_id(_matrix_indices character varying[], _bit character varying, output_param_id character varying) RETURNS character varying
    LANGUAGE plpgsql
    AS $$
DECLARE
        _sql varchar := 'SELECT pek1.param_id 
                        from param_piezo pek1, param_piezo pek2,
                        (select cast(%L  as varchar[]) as mi) x 
                        where pek1.param_id = pek2.param_id
                        and pek1.index_id < pek2.index_id
                        and cast(mi[array_position(mi,pek1.index_id)+1] as integer) = pek1.value
                        and cast(mi[array_position(mi,pek2.index_id)+1] as integer) = pek2.value
                        and pek1.param_matrix = %L';
        _m_i varchar[] := _matrix_indices;
        _out varchar;
BEGIN
        --- Flip the bit given as output bit
        IF _bit is not null THEN
        _m_i[array_position(_matrix_indices, _bit)+1] := cast(1-cast(_m_i[array_position(_matrix_indices, _bit)+1] as integer) as varchar);
        EXECUTE format(_sql, _m_i, output_param_id) INTO _out;
        return _out;
        ELSE
        RETURN NULL;
        END IF;
END
$$;


ALTER FUNCTION public.get_matrix_output_param_id(_matrix_indices character varying[], _bit character varying, output_param_id character varying) OWNER TO mena;

--
-- Name: integral_function(character varying, double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.integral_function(_function character varying, _params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql
    AS $$
declare  
        _out double precision[];
        _sql text := 'SELECT cast(%I(cast(%L as double precision[])) as double precision[]) as question';
Begin   
        
        EXECUTE format(_sql, _function, _params ) INTO _out;
        RETURN _out;
END
$$;


ALTER FUNCTION public.integral_function(_function character varying, _params double precision[]) OWNER TO mena;

--
-- Name: interpol_function(character varying, public.interpol_type[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.interpol_function(_function character varying, _params public.interpol_type[]) RETURNS real
    LANGUAGE plpgsql STABLE
    AS $$
declare  
	y real;
	x real;
	_out real;
	_sql text := 'SELECT cast(%I(cast(%L as interpol_type[])) as real) as question';
Begin   
	
	EXECUTE format(_sql, _function, _params ) INTO _out;
	RETURN _out;
END
$$;


ALTER FUNCTION public.interpol_function(_function character varying, _params public.interpol_type[]) OWNER TO mena;

--
-- Name: is_distinct(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.is_distinct(_params double precision[]) RETURNS boolean
    LANGUAGE plpgsql STABLE
    AS $$
DECLARE
out double precision;
sql text := 'Select count(*) from (select distinct * from unnest(cast(%L as double precision[]))) X';
BEGIN
EXECUTE format(sql,_params) into out;
Return cardinality(_params)=out;
END;
$$;


ALTER FUNCTION public.is_distinct(_params double precision[]) OWNER TO mena;

--
-- Name: is_length_1(character varying[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.is_length_1(_init character varying[]) RETURNS boolean
    LANGUAGE plpgsql STABLE
    AS $$
-- declare  
Begin   
RETURN array_length(_init,1)=1;
END
$$;


ALTER FUNCTION public.is_length_1(_init character varying[]) OWNER TO mena;

--
-- Name: linear_interpol_initial_stress(public.interpol_type[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.linear_interpol_initial_stress(_params public.interpol_type[]) RETURNS real
    LANGUAGE plpgsql
    AS $$
declare  
	out real;
	n int;
	x_g real;
	xy1_tmp interpol_type;
	xy2_tmp interpol_type;
	tmp interpol_type[];
	xy1 real[];
	xy2 real[];
	y_g real;
	row_content interpol_type;
Begin 
    
	SELECT x from unnest(_params) where y is null into x_g;
   

	SELECT x,y from unnest(_params) where y is not null and x = (Select max(x) from unnest(_params) where x <= x_g)  into xy1_tmp; 
    SELECT x,y from unnest(_params) where y is not null and x = (Select min(x) from unnest(_params) where x >= x_g)  into xy2_tmp; 
    
    
	if xy1_tmp is null or xy2_tmp is null THEN
		FOR row_content in
			SELECT x,y from unnest(_params) where y is not null order by abs(x-x_g) asc limit 2
		LOOP
			tmp := array_append(tmp, cast(row_content as interpol_type));
		END LOOP;
		xy1_tmp = tmp[2];
		xy2_tmp = tmp[3];
	END IF;
	
	xy1 := ARRAY[xy1_tmp.x, xy1_tmp.y];
	xy2 := ARRAY[xy2_tmp.x, xy2_tmp.y];
	

	if x_g=xy1[1] then
		y_g = xy1[2];
	elsif x_g=xy2[1] then
		y_g = xy2[2];
	else
		y_g = (xy1[2]-xy2[2])/(xy1[1]-xy2[1])*(x_g-xy1[1])+xy1[2];
	end if;
	
	RETURN y_g;
END
$$;


ALTER FUNCTION public.linear_interpol_initial_stress(_params public.interpol_type[]) OWNER TO mena;

--
-- Name: matrix_filter(public.matrix_filter_agg_type[], public.matrix_filter_agg_type); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.matrix_filter(oldrec public.matrix_filter_agg_type[], newrec public.matrix_filter_agg_type) RETURNS public.matrix_filter_agg_type[]
    LANGUAGE plpgsql STABLE
    AS $$
declare  
i integer;
Begin   

FOR i IN 1..COALESCE(array_length(oldrec,1),1)
LOOP
	IF (oldrec[i]).matr = (newrec).matr
	THEN
		IF (oldrec[i]).step > (newrec).step THEN
			oldrec[i] = newrec;
		END IF;
		RETURN oldrec;
	END IF;
END LOOP;
RETURN oldrec || newrec;
END
$$;


ALTER FUNCTION public.matrix_filter(oldrec public.matrix_filter_agg_type[], newrec public.matrix_filter_agg_type) OWNER TO mena;

--
-- Name: matrix_filter_end(public.matrix_filter_agg_type[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.matrix_filter_end(oldrec public.matrix_filter_agg_type[]) RETURNS integer
    LANGUAGE plpgsql STABLE
    AS $$
declare  
i integer;
Begin   

RETURN (SELECT max(step) FROM unnest(oldrec));
END
$$;


ALTER FUNCTION public.matrix_filter_end(oldrec public.matrix_filter_agg_type[]) OWNER TO mena;

--
-- Name: matrix_func_1(character varying[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.matrix_func_1(ar character varying[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$

import sympy as sy

def matrix_einlesen(ar):
    import numpy as np
    import sympy as sy

    var_e = int(ar[0][ar[0].index('e') + 1])
    var_m = int(ar[0][ar[0].index('m') + 1])

    A11 = None
    A12 = None
    x = sy.S('x')

    matrix_A = sy.Matrix([[x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [0., 0., 0., x, 0., 0.],
                          [0., 0., 0., 0., x, 0.], [0., 0., 0., 0., 0., x]])
    matrix_B = sy.Matrix([[0., 0., 0., 0., x, 0.], [0., 0., 0., x, 0., 0.], [x, x, x, 0., 0., 0.]])
    matrix_C = sy.Matrix([[x, 0., 0.], [0., x, 0.], [0., 0., x]])

    for i in range(0, len(ar)):
        ar_row = ar[i]
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_A[0, 0] = float(ar_row[5])
            matrix_A[1, 1] = float(ar_row[5])
            A11 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '2':
            matrix_A[0, 1] = float(ar_row[5])
            matrix_A[1, 0] = float(ar_row[5])
            A12 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '3':
            matrix_A[0, 2] = float(ar_row[5])
            matrix_A[2, 0] = float(ar_row[5])
            matrix_A[1, 2] = float(ar_row[5])
            matrix_A[2, 1] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_A[2, 2] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '5' and ar_row[4] == '5':
            matrix_A[4, 4] = float(ar_row[5])
            matrix_A[3, 3] = float(ar_row[5])
        if A11 is not None and A12 is not None:
            matrix_A[5, 5] = (A11 - A12) / (2 ** calc_vorzeichen(var_m))
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '1':
            matrix_B[2, 0] = float(ar_row[5])
            matrix_B[2, 1] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_B[2, 2] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '1' and ar_row[4] == '5':
            matrix_B[0, 4] = float(ar_row[5])
            matrix_B[1, 3] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_C[0, 0] = float(ar_row[5])
            matrix_C[1, 1] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_C[2, 2] = float(ar_row[5])

    return matrix_A, matrix_B, matrix_C

def matrix_rausschreiben(ar):
    import sympy as sy
    return_matrix = []
    if (ar.shape == (6, 6)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[0, 1] != 0. and type(ar[0, 1]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 2, ar[0, 1]])
        if (ar[0, 2] != 0. and type(ar[0, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 3, ar[0, 2]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[4, 4] != 0. and type(ar[4, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 5, 2, 5, ar[4, 4]])

    if (ar.shape == (3, 6)):
        if (ar[2, 0] != 0. and type(ar[2, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 1, ar[2, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[0, 4] != 0. and type(ar[0, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 5, ar[0, 4]])

    if (ar.shape == (3, 3)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])

    return return_matrix

def calc_vorzeichen(i):
    if i == 0:
        return -1
    elif i == 1:
        return 1
    else:
        return None

matrix_A1, matrix_B1, matrix_C1 = matrix_einlesen(ar)

x = sy.S('x')

if matrix_A1[0,0]==x or matrix_A1[0,1]==x or matrix_A1[0,2]==x or matrix_A1[2,2]==x:
	matrix_A = x*sy.ones(6,6)
	matrix_A[4,4] = 1/matrix_A1[4,4]
	matrix_A[3,3] = 1/matrix_A1[3,3]
	return matrix_rausschreiben(matrix_A)

if matrix_A1.det() == 0:
	return []
matrix_A = matrix_A1.inv()

return matrix_rausschreiben(matrix_A)
 $$;


ALTER FUNCTION public.matrix_func_1(ar character varying[]) OWNER TO mena;

--
-- Name: matrix_func_2(character varying[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.matrix_func_2(ar character varying[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$
import sympy as sy                                        
def matrix_einlesen(ar):
    import numpy as np
    import sympy as sy

    var_e = int(ar[0][ar[0].index('e') + 1])
    var_m = int(ar[0][ar[0].index('m') + 1])

    A11 = None
    A12 = None
    x = sy.S('x')

    matrix_A = sy.Matrix([[x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [0., 0., 0., x, 0., 0.],
                          [0., 0., 0., 0., x, 0.], [0., 0., 0., 0., 0., x]])
    matrix_B = sy.Matrix([[0., 0., 0., 0., x, 0.], [0., 0., 0., x, 0., 0.], [x, x, x, 0., 0., 0.]])
    matrix_C = sy.Matrix([[x, 0., 0.], [0., x, 0.], [0., 0., x]])

    for i in range(0, len(ar)):
        ar_row = ar[i]
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_A[0, 0] = float(ar_row[5])
            matrix_A[1, 1] = float(ar_row[5])
            A11 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '2':
            matrix_A[0, 1] = float(ar_row[5])
            matrix_A[1, 0] = float(ar_row[5])
            A12 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '3':
            matrix_A[0, 2] = float(ar_row[5])
            matrix_A[2, 0] = float(ar_row[5])
            matrix_A[1, 2] = float(ar_row[5])
            matrix_A[2, 1] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_A[2, 2] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '5' and ar_row[4] == '5':
            matrix_A[4, 4] = float(ar_row[5])
            matrix_A[3, 3] = float(ar_row[5])
        if A11 is not None and A12 is not None:
            matrix_A[5, 5] = (A11 - A12) / (2 ** calc_vorzeichen(var_m))
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '1':
            matrix_B[2, 0] = float(ar_row[5])
            matrix_B[2, 1] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_B[2, 2] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '1' and ar_row[4] == '5':
            matrix_B[0, 4] = float(ar_row[5])
            matrix_B[1, 3] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_C[0, 0] = float(ar_row[5])
            matrix_C[1, 1] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_C[2, 2] = float(ar_row[5])

    return matrix_A, matrix_B, matrix_C

def matrix_rausschreiben(ar):
    import sympy as sy
    return_matrix = []
    if (ar.shape == (6, 6)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[0, 1] != 0. and type(ar[0, 1]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 2, ar[0, 1]])
        if (ar[0, 2] != 0. and type(ar[0, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 3, ar[0, 2]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[4, 4] != 0. and type(ar[4, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 5, 2, 5, ar[4, 4]])

    if (ar.shape == (3, 6)):
        if (ar[2, 0] != 0. and type(ar[2, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 1, ar[2, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[0, 4] != 0. and type(ar[0, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 5, ar[0, 4]])

    if (ar.shape == (3, 3)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])

    return return_matrix

def calc_vorzeichen(i):
    if i == 0:
        return -1
    elif i == 1:
        return 1
    else:
        return None                         

matrix_A1, matrix_B1, matrix_C1 =matrix_einlesen(ar)               

x = sy.S('x')

if matrix_A1[0,0]==x or matrix_A1[0,1]==x or matrix_A1[0,2]==x or matrix_A1[2,2]==x:
     matrix_B = x*sy.ones(3,6)
     matrix_B[0,4] = matrix_B1[0,4]/matrix_A1[4,4]
     return matrix_rausschreiben(matrix_B)

                                                                        
matrix_B = matrix_B1 * (matrix_A1.inv())                               
                                                                      
return matrix_rausschreiben(matrix_B)$$;


ALTER FUNCTION public.matrix_func_2(ar character varying[]) OWNER TO mena;

--
-- Name: matrix_func_3(character varying[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.matrix_func_3(ar character varying[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$
import sympy as sy
import numpy as np

def matrix_einlesen(ar):
    import numpy as np
    import sympy as sy

    var_e = int(ar[0][ar[0].index('e') + 1])
    var_m = int(ar[0][ar[0].index('m') + 1])

    A11 = None
    A12 = None
    x = sy.S('x')

    matrix_A = sy.Matrix([[x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [0., 0., 0., x, 0., 0.],
                          [0., 0., 0., 0., x, 0.], [0., 0., 0., 0., 0., x]])
    matrix_B = sy.Matrix([[0., 0., 0., 0., x, 0.], [0., 0., 0., x, 0., 0.], [x, x, x, 0., 0., 0.]])
    matrix_C = sy.Matrix([[x, 0., 0.], [0., x, 0.], [0., 0., x]])

    for i in range(0, len(ar)):
        ar_row = ar[i]
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_A[0, 0] = float(ar_row[5])
            matrix_A[1, 1] = float(ar_row[5])
            A11 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '2':
            matrix_A[0, 1] = float(ar_row[5])
            matrix_A[1, 0] = float(ar_row[5])
            A12 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '3':
            matrix_A[0, 2] = float(ar_row[5])
            matrix_A[2, 0] = float(ar_row[5])
            matrix_A[1, 2] = float(ar_row[5])
            matrix_A[2, 1] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_A[2, 2] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '5' and ar_row[4] == '5':
            matrix_A[4, 4] = float(ar_row[5])
            matrix_A[3, 3] = float(ar_row[5])
        if A11 is not None and A12 is not None:
            matrix_A[5, 5] = (A11 - A12) / (2 ** calc_vorzeichen(var_m))
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '1':
            matrix_B[2, 0] = float(ar_row[5])
            matrix_B[2, 1] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_B[2, 2] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '1' and ar_row[4] == '5':
            matrix_B[0, 4] = float(ar_row[5])
            matrix_B[1, 3] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_C[0, 0] = float(ar_row[5])
            matrix_C[1, 1] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_C[2, 2] = float(ar_row[5])

    return matrix_A, matrix_B, matrix_C

def matrix_rausschreiben(ar):
    import sympy as sy
    return_matrix = []
    if (ar.shape == (6, 6)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[0, 1] != 0. and type(ar[0, 1]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 2, ar[0, 1]])
        if (ar[0, 2] != 0. and type(ar[0, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 3, ar[0, 2]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[4, 4] != 0. and type(ar[4, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 5, 2, 5, ar[4, 4]])

    if (ar.shape == (3, 6)):
        if (ar[2, 0] != 0. and type(ar[2, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 1, ar[2, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[0, 4] != 0. and type(ar[0, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 5, ar[0, 4]])

    if (ar.shape == (3, 3)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])

    return return_matrix

def calc_vorzeichen(i):
    if i == 0:
        return -1
    elif i == 1:
        return 1
    else:
        return None

matrix_A1, matrix_B1, matrix_C1 = matrix_einlesen(ar)
var_e = int(ar[0][ar[0].index('e')+1])
var_m = int(ar[0][ar[0].index('m')+1])

x = sy.S('x')

if matrix_A1[0,0]==x or matrix_A1[0,1]==x or matrix_A1[0,2]==x or matrix_A1[2,2]==x:
    matrix_A = sy.zeros(6,6)
    matrix_A[0:3,0:3] = x*sy.ones(3,3)
    matrix_A[5,5] = x
    matrix_A[4,4] = 1/matrix_A1[4,4]
    matrix_A[3,3] = 1/matrix_A1[3,3]
else:
    matrix_A = matrix_A1.inv()

matrix_C = matrix_C1 - calc_vorzeichen(var_e)*calc_vorzeichen(var_m)*matrix_B1*matrix_A*np.transpose(matrix_B1)

return matrix_rausschreiben(matrix_C)$$;


ALTER FUNCTION public.matrix_func_3(ar character varying[]) OWNER TO mena;

--
-- Name: matrix_func_4(character varying[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.matrix_func_4(ar character varying[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$import sympy as sy
import numpy as np

def matrix_einlesen(ar):
    var_e = int(ar[0][ar[0].index('e') + 1])
    var_m = int(ar[0][ar[0].index('m') + 1])

    A11 = None
    A12 = None
    x = sy.S('x')

    matrix_A = sy.Matrix([[x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [0., 0., 0., x, 0., 0.],
                          [0., 0., 0., 0., x, 0.], [0., 0., 0., 0., 0., x]])
    matrix_B = sy.Matrix([[0., 0., 0., 0., x, 0.], [0., 0., 0., x, 0., 0.], [x, x, x, 0., 0., 0.]])
    matrix_C = sy.Matrix([[x, 0., 0.], [0., x, 0.], [0., 0., x]])

    for i in range(0, len(ar)):
        ar_row = ar[i]
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_A[0, 0] = float(ar_row[5])
            matrix_A[1, 1] = float(ar_row[5])
            A11 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '2':
            matrix_A[0, 1] = float(ar_row[5])
            matrix_A[1, 0] = float(ar_row[5])
            A12 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '3':
            matrix_A[0, 2] = float(ar_row[5])
            matrix_A[2, 0] = float(ar_row[5])
            matrix_A[1, 2] = float(ar_row[5])
            matrix_A[2, 1] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_A[2, 2] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '5' and ar_row[4] == '5':
            matrix_A[4, 4] = float(ar_row[5])
            matrix_A[3, 3] = float(ar_row[5])
        if A11 is not None and A12 is not None:
            matrix_A[5, 5] = (A11 - A12) / (2 ** calc_vorzeichen(var_m))
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '1':
            matrix_B[2, 0] = float(ar_row[5])
            matrix_B[2, 1] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_B[2, 2] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '1' and ar_row[4] == '5':
            matrix_B[0, 4] = float(ar_row[5])
            matrix_B[1, 3] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_C[0, 0] = float(ar_row[5])
            matrix_C[1, 1] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_C[2, 2] = float(ar_row[5])

    return matrix_A, matrix_B, matrix_C

def matrix_rausschreiben(ar):
    return_matrix = []
    if (ar.shape == (6, 6)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[0, 1] != 0. and type(ar[0, 1]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 2, ar[0, 1]])
        if (ar[0, 2] != 0. and type(ar[0, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 3, ar[0, 2]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[4, 4] != 0. and type(ar[4, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 5, 2, 5, ar[4, 4]])

    if (ar.shape == (3, 6)):
        if (ar[2, 0] != 0. and type(ar[2, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 1, ar[2, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[0, 4] != 0. and type(ar[0, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 5, ar[0, 4]])

    if (ar.shape == (3, 3)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])

    return return_matrix

def calc_vorzeichen(i):
    if i == 0:
        return -1
    elif i == 1:
        return 1
    else:
        return None

matrix_A1, matrix_B1, matrix_C1 = matrix_einlesen(ar)
var_e = int(ar[0][ar[0].index('e')+1])
var_m = int(ar[0][ar[0].index('m')+1])
x = sy.S('x')
if matrix_C1.det() == 0 or matrix_C1.det() == x**3:
    return []
matrix_A = matrix_A1 - calc_vorzeichen(var_e)*calc_vorzeichen(var_m)*np.transpose(matrix_B1)*matrix_C1.inv()*matrix_B1

return matrix_rausschreiben(matrix_A)$$;


ALTER FUNCTION public.matrix_func_4(ar character varying[]) OWNER TO mena;

--
-- Name: matrix_func_5(character varying[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.matrix_func_5(ar character varying[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$import sympy as sy
    

def matrix_einlesen(ar):
    var_e = int(ar[0][ar[0].index('e') + 1])
    var_m = int(ar[0][ar[0].index('m') + 1])

    A11 = None
    A12 = None
    x = sy.S('x')

    matrix_A = sy.Matrix([[x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [0., 0., 0., x, 0., 0.],
                          [0., 0., 0., 0., x, 0.], [0., 0., 0., 0., 0., x]])
    matrix_B = sy.Matrix([[0., 0., 0., 0., x, 0.], [0., 0., 0., x, 0., 0.], [x, x, x, 0., 0., 0.]])
    matrix_C = sy.Matrix([[x, 0., 0.], [0., x, 0.], [0., 0., x]])

    for i in range(0, len(ar)):
        ar_row = ar[i]
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_A[0, 0] = float(ar_row[5])
            matrix_A[1, 1] = float(ar_row[5])
            A11 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '2':
            matrix_A[0, 1] = float(ar_row[5])
            matrix_A[1, 0] = float(ar_row[5])
            A12 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '3':
            matrix_A[0, 2] = float(ar_row[5])
            matrix_A[2, 0] = float(ar_row[5])
            matrix_A[1, 2] = float(ar_row[5])
            matrix_A[2, 1] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_A[2, 2] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '5' and ar_row[4] == '5':
            matrix_A[4, 4] = float(ar_row[5])
            matrix_A[3, 3] = float(ar_row[5])
        if A11 is not None and A12 is not None:
            matrix_A[5, 5] = (A11 - A12) / (2 ** calc_vorzeichen(var_m))
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '1':
            matrix_B[2, 0] = float(ar_row[5])
            matrix_B[2, 1] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_B[2, 2] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '1' and ar_row[4] == '5':
            matrix_B[0, 4] = float(ar_row[5])
            matrix_B[1, 3] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_C[0, 0] = float(ar_row[5])
            matrix_C[1, 1] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_C[2, 2] = float(ar_row[5])

    return matrix_A, matrix_B, matrix_C

def matrix_rausschreiben(ar):
    return_matrix = []
    if (ar.shape == (6, 6)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[0, 1] != 0. and type(ar[0, 1]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 2, ar[0, 1]])
        if (ar[0, 2] != 0. and type(ar[0, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 3, ar[0, 2]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[4, 4] != 0. and type(ar[4, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 5, 2, 5, ar[4, 4]])

    if (ar.shape == (3, 6)):
        if (ar[2, 0] != 0. and type(ar[2, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 1, ar[2, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[0, 4] != 0. and type(ar[0, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 5, ar[0, 4]])

    if (ar.shape == (3, 3)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])

    return return_matrix

def calc_vorzeichen(i):
    if i == 0:
        return -1
    elif i == 1:
        return 1
    else:
        return None    

matrix_A1, matrix_B1, matrix_C1 = matrix_einlesen(ar)
x = sy.S('x')
if matrix_C1.det() == 0 or matrix_C1.det() == x**3:
    return []
matrix_B =(matrix_C1.inv()) * matrix_B1 

return matrix_rausschreiben(matrix_B)$$;


ALTER FUNCTION public.matrix_func_5(ar character varying[]) OWNER TO mena;

--
-- Name: matrix_func_6(character varying[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.matrix_func_6(ar character varying[]) RETURNS double precision[]
    LANGUAGE plpython3u
    AS $$import sympy as sy
    

def matrix_einlesen(ar):
    var_e = int(ar[0][ar[0].index('e') + 1])
    var_m = int(ar[0][ar[0].index('m') + 1])

    A11 = None
    A12 = None
    x = sy.S('x')

    matrix_A = sy.Matrix([[x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [x, x, x, 0., 0., 0.], [0., 0., 0., x, 0., 0.],
                          [0., 0., 0., 0., x, 0.], [0., 0., 0., 0., 0., x]])
    matrix_B = sy.Matrix([[0., 0., 0., 0., x, 0.], [0., 0., 0., x, 0., 0.], [x, x, x, 0., 0., 0.]])
    matrix_C = sy.Matrix([[x, 0., 0.], [0., x, 0.], [0., 0., x]])

    for i in range(0, len(ar)):
        ar_row = ar[i]
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_A[0, 0] = float(ar_row[5])
            matrix_A[1, 1] = float(ar_row[5])
            A11 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '2':
            matrix_A[0, 1] = float(ar_row[5])
            matrix_A[1, 0] = float(ar_row[5])
            A12 = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '1' and ar_row[4] == '3':
            matrix_A[0, 2] = float(ar_row[5])
            matrix_A[2, 0] = float(ar_row[5])
            matrix_A[1, 2] = float(ar_row[5])
            matrix_A[2, 1] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_A[2, 2] = float(ar_row[5])
        if ar_row[0] == 'A' and ar_row[2] == '5' and ar_row[4] == '5':
            matrix_A[4, 4] = float(ar_row[5])
            matrix_A[3, 3] = float(ar_row[5])
        if A11 is not None and A12 is not None:
            matrix_A[5, 5] = (A11 - A12) / (2 ** calc_vorzeichen(var_m))
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '1':
            matrix_B[2, 0] = float(ar_row[5])
            matrix_B[2, 1] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_B[2, 2] = float(ar_row[5])
        if ar_row[0] == 'B' and ar_row[2] == '1' and ar_row[4] == '5':
            matrix_B[0, 4] = float(ar_row[5])
            matrix_B[1, 3] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '1' and ar_row[4] == '1':
            matrix_C[0, 0] = float(ar_row[5])
            matrix_C[1, 1] = float(ar_row[5])
        if ar_row[0] == 'C' and ar_row[2] == '3' and ar_row[4] == '3':
            matrix_C[2, 2] = float(ar_row[5])

    return matrix_A, matrix_B, matrix_C

def matrix_rausschreiben(ar):
    return_matrix = []
    if (ar.shape == (6, 6)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[0, 1] != 0. and type(ar[0, 1]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 2, ar[0, 1]])
        if (ar[0, 2] != 0. and type(ar[0, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 3, ar[0, 2]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[4, 4] != 0. and type(ar[4, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 5, 2, 5, ar[4, 4]])

    if (ar.shape == (3, 6)):
        if (ar[2, 0] != 0. and type(ar[2, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 1, ar[2, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])
        if (ar[0, 4] != 0. and type(ar[0, 4]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 5, ar[0, 4]])

    if (ar.shape == (3, 3)):
        if (ar[0, 0] != 0. and type(ar[0, 0]) == sy.core.numbers.Float):
            return_matrix.append([1, 1, 2, 1, ar[0, 0]])
        if (ar[2, 2] != 0. and type(ar[2, 2]) == sy.core.numbers.Float):
            return_matrix.append([1, 3, 2, 3, ar[2, 2]])

    return return_matrix

def calc_vorzeichen(i):
    if i == 0:
        return -1
    elif i == 1:
        return 1
    else:
        return None    

matrix_A1, matrix_B1, matrix_C1 = matrix_einlesen(ar)
x = sy.S('x')
if matrix_C1.det() == 0 or matrix_C1.det() == x**3:
    return []
matrix_C = matrix_C1.inv()

return matrix_rausschreiben(matrix_C)$$;


ALTER FUNCTION public.matrix_func_6(ar character varying[]) OWNER TO mena;

--
-- Name: matrix_function(character varying, character varying[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.matrix_function(_function character varying, _params character varying[]) RETURNS double precision[]
    LANGUAGE plpgsql STABLE
    AS $$
declare  
	_out double precision[];
	_sql text := 'SELECT cast(%I(cast(%L as varchar[])) as double precision[]) as question';
Begin   
	
	EXECUTE format(_sql, _function, _params ) INTO _out;
	RETURN _out;
END
$$;


ALTER FUNCTION public.matrix_function(_function character varying, _params character varying[]) OWNER TO mena;

--
-- Name: mech_mod_v1(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.mech_mod_v1(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql STABLE
    AS $$                                              
 declare                                          
  eps double precision = _params[1];
  E_A double precision = _params[2];
  E_M double precision = _params[3];
 xi double precision = _params[4];
 Begin 
       RETURN ARRAY[eps * (E_A - (E_A-E_M)*xi)];
 
 END                                             
 $$;


ALTER FUNCTION public.mech_mod_v1(_params double precision[]) OWNER TO mena;

--
-- Name: mech_mod_v2(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.mech_mod_v2(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql STABLE
    AS $$                                              
 declare                                          
  eps double precision = _params[1];
  E_A double precision = _params[2];
  E_M double precision = _params[3];
 xi double precision = _params[4];
 eps_ps double precision = _params[5];
 H double precision = _params[6];
 Begin 
       RETURN ARRAY[eps_ps * (E_A - (E_A-E_M)*xi)+ (eps-eps_ps)*H];
 
 END                                             
 $$;


ALTER FUNCTION public.mech_mod_v2(_params double precision[]) OWNER TO mena;

--
-- Name: mech_mod_v3(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.mech_mod_v3(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql STABLE
    AS $$                                              
 declare                                          
  eps double precision = _params[1];
  E_A double precision = _params[2];
  E_M double precision = _params[3];
 xi double precision = _params[4];
 eps_ps double precision = _params[5];
 H double precision = _params[6];
 eps_pf double precision = _params[7];
 E_D double precision = _params[8];  
 Begin 
       RETURN ARRAY[eps_ps * (E_A - (E_A-E_M)*xi)+ (eps-eps_ps)*H + (eps - eps_pf)*E_D];
 
 END                                             
 $$;


ALTER FUNCTION public.mech_mod_v3(_params double precision[]) OWNER TO mena;

--
-- Name: mod_el_druck(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.mod_el_druck(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql IMMUTABLE
    AS $$                                                             
  declare                                                                  
   rel_perm double precision = _params[1];                                 
   vac_perm double precision = _params[2];                            
   el_field double precision = _params[3];                                  
                                                                           
  Begin                                                                    
        RETURN ARRAY[rel_perm * vac_perm * pow(el_field,2)];                         
                                                                           
  end                                                                      
  $$;


ALTER FUNCTION public.mod_el_druck(_params double precision[]) OWNER TO mena;

--
-- Name: mod_el_pres(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.mod_el_pres(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql IMMUTABLE
    AS $$                                                             
  declare                                                                  
   rel_perm double precision = _params[1];                                 
   vac_perm double precision = _params[2];                            
   el_field double precision = _params[3];                                  
                                                                           
  Begin                                                                    
        RETURN ARRAY[rel_perm * vac_perm * pow(el_field,2)];                         
                                                                           
  end                                                                      
  $$;


ALTER FUNCTION public.mod_el_pres(_params double precision[]) OWNER TO mena;

--
-- Name: round9(double precision); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.round9(input double precision) RETURNS double precision
    LANGUAGE plpgsql STABLE
    AS $$                                                   
         DECLARE                                                 
                                                                 
         BEGIN                                                   
          IF abs(input) >= 100 THEN                                   
                 RETURN round(input::decimal,0);                 
          ELSIF abs(input) > 10 THEN                                  
                 RETURN round(input::decimal,1);                 
          ELSIF abs(input) >= 1 THEN                                  
                 RETURN round(input::decimal,2);                 
          ELSIF abs(input) > 0.1 THEN                                 
                 RETURN round(input::decimal,3);                 
           ELSIF abs(input) >= 0.01 THEN                              
                 RETURN round(input::decimal,4);                 
          ELSE                                                   
                 RETURN input;                                   
         END IF;                                                 
         END                                                     
 $$;


ALTER FUNCTION public.round9(input double precision) OWNER TO mena;

--
-- Name: simple_function(character varying, double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.simple_function(_function character varying, _params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql STABLE
    AS $$                                                                                                  
  declare                                                                                                       
          _out double precision [];                                                                                
          _sql text := 'SELECT cast(%I(cast(%L as double precision[])) as double precision[]) as question';     
  Begin                                                                                                         
                                                                                                                
          EXECUTE format(_sql, _function, _params ) INTO _out;                                                  
          RETURN _out;                                                                                          
  END                                                                                                           
  $$;


ALTER FUNCTION public.simple_function(_function character varying, _params double precision[]) OWNER TO mena;

--
-- Name: size_1_or_2nd_diff(character varying[], character varying[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.size_1_or_2nd_diff(_init character varying[], _ar character varying[]) RETURNS character varying[]
    LANGUAGE plpgsql STABLE
    AS $$
-- declare  
Begin
IF array_length(_ar,1)<1 OR array_length(_init,1)<1 THEN
	RETURN ARRAY[];
ELSIF array_length(_ar,1)=1 OR array_length(_init,1)=1 OR NOT _ar[2] = _init[2] THEN
	RETURN ARRAY['fertig'];
ELSE
	RETURN _ar;
END IF;
END
$$;


ALTER FUNCTION public.size_1_or_2nd_diff(_init character varying[], _ar character varying[]) OWNER TO mena;

--
-- Name: test_compare_condition(double precision[], public.operator_types[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.test_compare_condition(_params double precision[], _operators public.operator_types[]) RETURNS boolean
    LANGUAGE plpgsql
    AS $$ 

Begin 
if array_length(_params, 1) is null 
    then return false;
end if;

if MOD(array_length(_params, 1), 2) = 1 then 
    return false;
end if;

For i in 1..array_length(_params, 1) by 2 loop 
    If not test_condition(_params[i], _params[i+1], _operators[i]) then 
        return false;
    End if;
End loop;

Return true;

End;

$$;


ALTER FUNCTION public.test_compare_condition(_params double precision[], _operators public.operator_types[]) OWNER TO mena;

--
-- Name: test_condition(double precision, double precision, public.operator_types); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.test_condition(_value double precision, _condition double precision, _operator public.operator_types) RETURNS boolean
    LANGUAGE plpgsql
    AS $$
DECLARE
sql text; 
out boolean;
BEGIN
sql := 'Select cast(%L as double precision)'||_operator||'cast(%L as double precision)';
EXECUTE format(sql,_value,_condition) INTO out;
RETURN out;
END
$$;


ALTER FUNCTION public.test_condition(_value double precision, _condition double precision, _operator public.operator_types) OWNER TO mena;

--
-- Name: unnest_nd_1d(anyarray); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.unnest_nd_1d(a anyarray, OUT a_1d anyarray) RETURNS SETOF anyarray
    LANGUAGE plpgsql IMMUTABLE STRICT PARALLEL SAFE
    AS $$
BEGIN
   FOREACH a_1d SLICE 1 IN ARRAY a LOOP
      RETURN NEXT;
   END LOOP;
END
$$;


ALTER FUNCTION public.unnest_nd_1d(a anyarray, OUT a_1d anyarray) OWNER TO mena;

--
-- Name: ym_init_slope(double precision[]); Type: FUNCTION; Schema: public; Owner: mena
--

CREATE FUNCTION public.ym_init_slope(_params double precision[]) RETURNS double precision[]
    LANGUAGE plpgsql STABLE
    AS $$                                              
 declare                                          
  sigma real = _params[1];
  epsilon real = _params[2];
 Begin 
 If epsilon = 0 THEN
 return NULL;
 else
       RETURN array[sigma / epsilon];
       End if;
 END                                             
 $$;


ALTER FUNCTION public.ym_init_slope(_params double precision[]) OWNER TO mena;

--
-- Name: matrix_filter_agg(public.matrix_filter_agg_type); Type: AGGREGATE; Schema: public; Owner: mena
--

CREATE AGGREGATE public.matrix_filter_agg(mat public.matrix_filter_agg_type) (
    SFUNC = public.matrix_filter,
    STYPE = public.matrix_filter_agg_type[],
    INITCOND = '{}',
    FINALFUNC = public.matrix_filter_end
);


ALTER AGGREGATE public.matrix_filter_agg(mat public.matrix_filter_agg_type) OWNER TO mena;

--
-- Name: size_1_or_2nd_diff(character varying[]); Type: AGGREGATE; Schema: public; Owner: mena
--

CREATE AGGREGATE public.size_1_or_2nd_diff(params character varying[]) (
    SFUNC = public.size_1_or_2nd_diff,
    STYPE = character varying[],
    INITCOND = '{}',
    FINALFUNC = public.is_length_1
);


ALTER AGGREGATE public.size_1_or_2nd_diff(params character varying[]) OWNER TO mena;

SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: characteristic_curve; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.characteristic_curve (
    curve_id character varying NOT NULL,
    curve_name character varying,
    curve_descripton character varying
);


ALTER TABLE public.characteristic_curve OWNER TO mena;

--
-- Name: component; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.component (
    component_id character varying NOT NULL,
    component_name character varying
);


ALTER TABLE public.component OWNER TO mena;

--
-- Name: component_parts; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.component_parts (
    concr_component_id character varying NOT NULL,
    has_part character varying NOT NULL
);


ALTER TABLE public.component_parts OWNER TO mena;

--
-- Name: concr_component; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.concr_component (
    concr_component_id character varying NOT NULL,
    component_id character varying,
    concr_component_name character varying
);


ALTER TABLE public.concr_component OWNER TO mena;

--
-- Name: concr_curve_params; Type: VIEW; Schema: public; Owner: mena
--

CREATE VIEW public.concr_curve_params AS
 SELECT (public.function_concr_curve_params()).id_tuple AS id_tuple,
    (public.function_concr_curve_params()).id_curve AS id_curve,
    (public.function_concr_curve_params()).curve_id AS curve_id,
    (public.function_concr_curve_params()).mat_sample_id AS mat_sample_id,
    (public.function_concr_curve_params()).meas_time AS meas_time;


ALTER TABLE public.concr_curve_params OWNER TO mena;

--
-- Name: curve_params; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.curve_params (
    curve_id character varying NOT NULL,
    param_id character varying NOT NULL
);


ALTER TABLE public.curve_params OWNER TO mena;

--
-- Name: concr_curve; Type: VIEW; Schema: public; Owner: mena
--

CREATE VIEW public.concr_curve AS
 SELECT DISTINCT (curve_params.concr_curve_id || concr_curve_params.id_curve) AS concr_curve_id,
    curve_params.curve_id,
    concr_curve_params.id_curve,
    curve_params.param_id,
    curve_params.in_number
   FROM public.concr_curve_params,
    ( SELECT cp.concr_curve_id,
            cp.curve_id,
            cp.param1 AS param_id,
            1 AS in_number
           FROM ( SELECT (((c.curve_id)::text || (c.param_id)::text) || (d.param_id)::text) AS concr_curve_id,
                    c.curve_id,
                    c.param_id AS param1,
                    d.param_id AS param2
                   FROM public.curve_params c,
                    public.curve_params d
                  WHERE (((c.curve_id)::text = (d.curve_id)::text) AND ((c.param_id)::text <> (d.param_id)::text))) cp
        UNION
         SELECT cp.concr_curve_id,
            cp.curve_id,
            cp.param2 AS param_id,
            2 AS in_number
           FROM ( SELECT (((c.curve_id)::text || (c.param_id)::text) || (d.param_id)::text) AS concr_curve_id,
                    c.curve_id,
                    c.param_id AS param1,
                    d.param_id AS param2
                   FROM public.curve_params c,
                    public.curve_params d
                  WHERE (((c.curve_id)::text = (d.curve_id)::text) AND ((c.param_id)::text <> (d.param_id)::text))) cp) curve_params
  WHERE ((concr_curve_params.curve_id)::text = (curve_params.curve_id)::text);


ALTER TABLE public.concr_curve OWNER TO mena;

--
-- Name: concr_curve_params_temp; Type: VIEW; Schema: public; Owner: mena
--

CREATE VIEW public.concr_curve_params_temp AS
 SELECT (public.function_concr_curve_params_temp()).concr_tuple_id AS concr_tuple_id,
    (public.function_concr_curve_params_temp()).id_tuple AS id_tuple,
    (public.function_concr_curve_params_temp()).concr_curve_id AS concr_curve_id,
    (public.function_concr_curve_params_temp()).concr_param_id1 AS concr_param_id1,
    (public.function_concr_curve_params_temp()).concr_param_id2 AS concr_param_id2;


ALTER TABLE public.concr_curve_params_temp OWNER TO mena;

--
-- Name: concr_meas; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.concr_meas (
    concr_meas_id character varying NOT NULL,
    meas_id character varying NOT NULL,
    meas_date date,
    description character varying
);


ALTER TABLE public.concr_meas OWNER TO mena;

--
-- Name: concr_model; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.concr_model (
    concr_model_id character varying NOT NULL,
    mod_id character varying NOT NULL,
    description character varying
);


ALTER TABLE public.concr_model OWNER TO mena;

--
-- Name: concr_param; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.concr_param (
    concr_param_id character varying NOT NULL,
    mod_meas_id character varying,
    param_id character varying NOT NULL,
    mat_sample_id character varying NOT NULL,
    meas_time integer,
    value double precision
);


ALTER TABLE public.concr_param OWNER TO mena;

--
-- Name: concr_param_matrix; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.concr_param_matrix (
    concr_param_id character varying NOT NULL,
    index_id integer NOT NULL,
    value integer NOT NULL
);


ALTER TABLE public.concr_param_matrix OWNER TO mena;

--
-- Name: concr_param_matrix_additional; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.concr_param_matrix_additional (
    concr_param_id character varying NOT NULL,
    index_id integer NOT NULL,
    value integer,
    callversion integer
);


ALTER TABLE public.concr_param_matrix_additional OWNER TO mena;

--
-- Name: concr_param_view; Type: VIEW; Schema: public; Owner: mena
--

CREATE VIEW public.concr_param_view AS
 SELECT (public.concr_param_function(0, 20, 0)).concr_param_id AS concr_param_id,
    (public.concr_param_function(0, 20, 0)).mod_meas_id AS mod_meas_id,
    (public.concr_param_function(0, 20, 0)).param_id AS param_id,
    (public.concr_param_function(0, 20, 0)).mat_sample_id AS mat_sample_id,
    (public.concr_param_function(0, 20, 0)).meas_time AS meas_time,
    (public.concr_param_function(0, 20, 0)).value AS value,
    (public.concr_param_function(0, 20, 0)).mod_id AS mod_id,
    (public.concr_param_function(0, 20, 0)).step AS step;


ALTER TABLE public.concr_param_view OWNER TO mena;

--
-- Name: constants; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.constants (
    param_id character varying NOT NULL,
    value numeric NOT NULL
);


ALTER TABLE public.constants OWNER TO mena;

--
-- Name: concr_parameter; Type: MATERIALIZED VIEW; Schema: public; Owner: mena
--

CREATE MATERIALIZED VIEW public.concr_parameter AS
 SELECT concr_param.concr_param_id,
    concr_param.mod_meas_id,
    concr_param.param_id,
    concr_param.mat_sample_id,
    concr_param.meas_time,
    concr_param.value,
    concr_mod_meas.mod_id AS mod_meas_id_abstract,
    0 AS step
   FROM (public.concr_param
     LEFT JOIN ( SELECT concr_model.concr_model_id AS mod_meas_id,
            concr_model.mod_id
           FROM public.concr_model
        UNION
         SELECT concr_meas.concr_meas_id AS mod_meas_id,
            concr_meas.meas_id AS mod_id
           FROM public.concr_meas) concr_mod_meas ON (((concr_mod_meas.mod_meas_id)::text = (concr_param.mod_meas_id)::text)))
UNION
 SELECT concr_param_view.concr_param_id,
    concr_param_view.mod_meas_id,
    concr_param_view.param_id,
    concr_param_view.mat_sample_id,
    concr_param_view.meas_time,
    concr_param_view.value,
    concr_param_view.mod_id AS mod_meas_id_abstract,
    concr_param_view.step
   FROM public.concr_param_view
  WHERE (concr_param_view.step >= 1)
UNION
 SELECT constants.param_id AS concr_param_id,
    NULL::character varying AS mod_meas_id,
    constants.param_id,
    NULL::character varying AS mat_sample_id,
    NULL::integer AS meas_time,
    constants.value,
    NULL::character varying AS mod_meas_id_abstract,
    0 AS step
   FROM public.constants
  WITH NO DATA;


ALTER TABLE public.concr_parameter OWNER TO mena;

--
-- Name: curve_condition; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.curve_condition (
    curve_id character varying NOT NULL,
    param_id character varying NOT NULL,
    value numeric,
    in_number integer NOT NULL,
    operator public.operator_types
);


ALTER TABLE public.curve_condition OWNER TO mena;

--
-- Name: material; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.material (
    mat_id character varying NOT NULL,
    mat_name character varying,
    mat_type character varying
);


ALTER TABLE public.material OWNER TO mena;

--
-- Name: material_condition; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.material_condition (
    mod_id character varying NOT NULL,
    mat_type character varying NOT NULL,
    abstract_object_id character varying
);


ALTER TABLE public.material_condition OWNER TO mena;

--
-- Name: material_types; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.material_types (
    mat_type character varying NOT NULL
);


ALTER TABLE public.material_types OWNER TO mena;

--
-- Name: measinput; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.measinput (
    meas_id character varying NOT NULL,
    param_id character varying NOT NULL
);


ALTER TABLE public.measinput OWNER TO mena;

--
-- Name: measoutput; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.measoutput (
    meas_id character varying NOT NULL,
    param_id character varying NOT NULL
);


ALTER TABLE public.measoutput OWNER TO mena;

--
-- Name: measurement; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.measurement (
    meas_id character varying NOT NULL,
    meas_name character varying
);


ALTER TABLE public.measurement OWNER TO mena;

--
-- Name: model; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.model (
    mod_id character varying NOT NULL,
    mod_name character varying,
    mod_type character varying,
    mod_equation character varying,
    function_name character varying
);


ALTER TABLE public.model OWNER TO mena;

--
-- Name: model_condition; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.model_condition (
    mod_id character varying NOT NULL,
    param_id character varying NOT NULL,
    value numeric,
    in_number integer,
    operator public.operator_types
);


ALTER TABLE public.model_condition OWNER TO mena;

--
-- Name: model_condition_compare; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.model_condition_compare (
    mod_id character varying NOT NULL,
    param_id character varying NOT NULL,
    in_number integer NOT NULL,
    operator public.operator_types NOT NULL,
    compare_position integer NOT NULL
);


ALTER TABLE public.model_condition_compare OWNER TO mena;

--
-- Name: model_condition_interpol; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.model_condition_interpol (
    mod_id character varying,
    param_id character varying,
    value numeric,
    in_number integer
);


ALTER TABLE public.model_condition_interpol OWNER TO mena;

--
-- Name: modelinput; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.modelinput (
    mod_id character varying NOT NULL,
    param_id character varying NOT NULL,
    in_number integer
);


ALTER TABLE public.modelinput OWNER TO mena;

--
-- Name: modelinput_matrix; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.modelinput_matrix (
    mod_id character varying NOT NULL,
    param_id character varying NOT NULL,
    index_id integer NOT NULL,
    value integer
);


ALTER TABLE public.modelinput_matrix OWNER TO mena;

--
-- Name: modeloutput; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.modeloutput (
    mod_id character varying NOT NULL,
    param_id character varying NOT NULL,
    out_number integer DEFAULT 1 NOT NULL
);


ALTER TABLE public.modeloutput OWNER TO mena;

--
-- Name: modeloutput_pek_bits; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.modeloutput_pek_bits (
    mod_id character varying NOT NULL,
    "bit" character varying
);


ALTER TABLE public.modeloutput_pek_bits OWNER TO mena;

--
-- Name: modeltype; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.modeltype (
    modeltype character varying NOT NULL
);


ALTER TABLE public.modeltype OWNER TO mena;

--
-- Name: param; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.param (
    param_id character varying NOT NULL,
    param_name character varying,
    param_unit character varying NOT NULL,
    param_symbol character varying NOT NULL
);


ALTER TABLE public.param OWNER TO mena;

--
-- Name: param_piezo; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.param_piezo (
    param_id character varying NOT NULL,
    index_id character(1) NOT NULL,
    param_matrix character varying,
    value double precision NOT NULL
);


ALTER TABLE public.param_piezo OWNER TO mena;

--
-- Name: param_semantic; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.param_semantic (
    param_id character varying NOT NULL,
    value double precision NOT NULL,
    semantic character varying
);


ALTER TABLE public.param_semantic OWNER TO mena;

--
-- Name: sample; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.sample (
    sample_id character varying NOT NULL,
    mat_id character varying NOT NULL,
    description character varying,
    sample_name character varying
);


ALTER TABLE public.sample OWNER TO mena;

--
-- Name: temp_step_view; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.temp_step_view (
    concr_param_id character varying,
    mod_meas_id character varying,
    param_id character varying,
    mat_sample_id character varying,
    meas_time integer,
    value double precision,
    mod_id character varying
);


ALTER TABLE public.temp_step_view OWNER TO mena;

--
-- Name: temp_view; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.temp_view (
    concr_param_id character varying,
    mod_meas_id character varying,
    param_id character varying,
    mat_sample_id character varying,
    meas_time integer,
    value double precision,
    mod_id character varying,
    step integer
);


ALTER TABLE public.temp_view OWNER TO mena;

--
-- Name: temp_view_results; Type: TABLE; Schema: public; Owner: mena
--

CREATE TABLE public.temp_view_results (
    concr_param_id character varying,
    mod_meas_id character varying,
    param_id character varying,
    mat_sample_id character varying,
    meas_time integer,
    value double precision,
    mod_id character varying,
    step integer
);


ALTER TABLE public.temp_view_results OWNER TO mena;

--
-- Data for Name: characteristic_curve; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.characteristic_curve (curve_id, curve_name, curve_descripton) FROM stdin;
magnetic_curve_easy_axis	magnetic curve easy axis	\N
magnetic_curve_hard_axis	magnetic curve hard axis	\N
hysteresis_martensite_content	hysterese model for martensit content over temperature	\N
stress-strain_curve_martensite	stress-strain curve in dependence of martensite	\N
stress-strain_curve_elastomer	stress-strain curve for elastomer materials	\N
curve_magnetic_stress	curve of magnetic stress	\N
\.


--
-- Data for Name: component; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.component (component_id, component_name) FROM stdin;
wire_actuator	wire actuator
stress-strain-stick	MSMA-stick
\.


--
-- Data for Name: component_parts; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.component_parts (concr_component_id, has_part) FROM stdin;
stick_Ni2MnGa_sample_1	Ni2MnGa_sample_1
\.


--
-- Data for Name: concr_component; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.concr_component (concr_component_id, component_id, concr_component_name) FROM stdin;
stick_Ni2MnGa_sample_1	stress-strain-stick	stick_Ni2MnGa_sample_1
wire_actuator_1	wire_actuator	wire_actuator_1
wire_actuator_NiTi#6	wire_actuator	wire_actuator_NiTi#6
\.


--
-- Data for Name: concr_meas; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.concr_meas (concr_meas_id, meas_id, meas_date, description) FROM stdin;
Ni2MnGa_sample1_meas_flux_over_field_1	meas_flux_over_field	2022-05-18	measurement of flux density over magnetic field, hard axis
Ni2MnGa_sample1_meas_flux_over_field_2	meas_flux_over_field	2022-05-18	measurement of flux density over magnetic field, easy axis
tens_1	tensile_test	2022-11-10	\N
elas1_measured_at	measured_at	\N	measured at Fraunhofer IAP
elas1_base_value	base_value	\N	base value for initial slope identification
dummy-meas	dummy-data	\N	dummy data
\.


--
-- Data for Name: concr_model; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.concr_model (concr_model_id, mod_id, description) FROM stdin;
concr_calc_max_blocking_stress_PC	calc_max_blocking_stress_PC	calculation of maxmimum blocking stress
\.


--
-- Data for Name: concr_param; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.concr_param (concr_param_id, mod_meas_id, param_id, mat_sample_id, meas_time, value) FROM stdin;
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_1	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	1	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_2	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	2	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_3	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	3	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_4	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	4	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_5	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	5	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_6	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	6	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_7	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	7	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_8	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	8	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_9	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	9	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_10	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	10	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_11	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	11	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_12	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	12	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_13	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	13	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_14	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	14	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_15	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	15	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_16	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	16	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_17	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	17	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_18	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	18	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_19	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	19	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_20	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	20	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_21	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	21	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_22	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	22	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_23	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	23	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_24	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	24	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_25	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	25	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_26	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	26	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_27	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	27	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_28	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	28	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_29	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	29	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_30	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	30	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_31	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	31	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_32	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	32	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_33	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	33	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_34	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	34	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_35	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	35	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_36	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	36	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_37	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	37	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_38	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	38	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_39	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	39	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_40	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	40	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_41	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	41	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_42	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	42	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_43	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	43	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_44	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	44	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_45	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	45	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_46	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	46	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_47	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	47	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_48	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	48	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_49	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	49	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_50	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	50	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_51	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	51	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_52	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	52	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_53	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	53	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_54	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	54	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_55	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	55	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_56	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	56	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_57	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	57	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_58	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	58	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_59	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	59	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_60	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	60	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_61	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	61	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_62	Ni2MnGa_sample1_meas_flux_over_field_1	easy_hard_axis	Ni2MnGa_sample_1	62	1
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_101	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	101	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_102	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	102	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_103	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	103	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_104	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	104	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_105	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	105	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_106	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	106	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_107	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	107	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_108	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	108	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_109	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	109	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_110	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	110	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_111	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	111	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_112	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	112	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_113	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	113	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_114	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	114	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_115	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	115	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_116	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	116	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_117	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	117	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_118	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	118	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_119	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	119	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_120	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	120	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_121	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	121	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_122	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	122	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_123	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	123	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_124	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	124	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_125	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	125	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_126	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	126	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_127	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	127	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_128	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	128	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_129	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	129	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_130	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	130	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_131	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	131	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_132	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	132	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_133	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	133	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_134	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	134	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_135	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	135	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_136	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	136	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_137	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	137	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_138	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	138	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_139	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	139	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_140	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	140	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_141	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	141	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_142	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	142	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_143	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	143	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_144	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	144	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_145	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	145	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_146	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	146	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_147	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	147	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_148	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	148	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_149	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	149	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_150	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	150	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_151	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	151	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_152	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	152	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_153	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	153	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_154	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	154	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_155	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	155	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_156	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	156	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_157	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	157	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_158	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	158	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_159	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	159	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_160	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	160	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_161	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	161	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_162	Ni2MnGa_sample1_meas_flux_over_field_2	easy_hard_axis	Ni2MnGa_sample_1	162	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_46	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	46	450000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_47	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	47	460000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_48	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	48	470000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_49	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	49	480000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_50	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	50	490000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_51	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	51	500000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_52	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	52	510000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_53	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	53	520000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_54	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	54	530000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_1	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	1	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_2	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	2	10000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_3	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	3	20000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_4	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	4	30000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_5	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	5	40000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_6	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	6	50000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_7	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	7	60000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_8	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	8	70000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_9	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	9	80000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_10	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	10	90000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_11	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	11	100000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_12	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	12	110000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_13	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	13	120000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_14	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	14	130000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_15	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	15	140000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_16	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	16	150000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_17	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	17	160000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_18	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	18	170000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_19	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	19	180000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_20	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	20	190000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_21	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	21	200000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_22	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	22	210000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_23	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	23	220000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_24	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	24	230000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_25	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	25	240000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_26	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	26	250000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_27	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	27	260000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_28	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	28	270000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_29	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	29	280000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_30	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	30	290000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_31	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	31	300000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_32	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	32	310000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_33	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	33	320000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_34	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	34	330000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_35	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	35	340000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_36	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	36	350000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_37	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	37	360000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_38	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	38	370000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_39	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	39	380000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_40	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	40	390000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_41	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	41	400000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_42	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	42	410000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_43	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	43	420000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_44	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	44	430000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_45	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	45	440000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_55	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	55	540000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_56	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	56	550000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_57	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	57	560000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_58	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	58	570000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_59	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	59	580000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_60	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	60	590000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_61	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	61	600000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_62	Ni2MnGa_sample1_meas_flux_over_field_1	magnetic_field_strength	Ni2MnGa_sample_1	62	610000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_101	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	101	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_102	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	102	10000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_103	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	103	20000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_104	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	104	30000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_105	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	105	40000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_106	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	106	50000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_107	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	107	60000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_108	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	108	70000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_109	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	109	80000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_110	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	110	90000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_111	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	111	100000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_112	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	112	110000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_113	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	113	120000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_114	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	114	130000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_115	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	115	140000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_116	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	116	150000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_117	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	117	160000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_118	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	118	170000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_119	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	119	180000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_120	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	120	190000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_121	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	121	200000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_122	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	122	210000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_123	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	123	220000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_124	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	124	230000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_125	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	125	240000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_126	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	126	250000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_127	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	127	260000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_128	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	128	270000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_129	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	129	280000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_130	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	130	290000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_131	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	131	300000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_132	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	132	310000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_133	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	133	320000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_134	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	134	330000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_135	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	135	340000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_136	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	136	350000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_137	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	137	360000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_138	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	138	370000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_139	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	139	380000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_140	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	140	390000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_141	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	141	400000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_142	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	142	410000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_143	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	143	420000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_144	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	144	430000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_145	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	145	440000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_146	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	146	450000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_147	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	147	460000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_148	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	148	470000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_149	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	149	480000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_150	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	150	490000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_151	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	151	500000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_152	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	152	510000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_153	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	153	520000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_154	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	154	530000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_155	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	155	540000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_156	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	156	550000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_157	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	157	560000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_158	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	158	570000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_159	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	159	580000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_160	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	160	590000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_161	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	161	600000
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_162	Ni2MnGa_sample1_meas_flux_over_field_2	magnetic_field_strength	Ni2MnGa_sample_1	162	610000
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_140	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	140	1.1593771
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_101	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	101	0
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_102	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	102	0.4355742
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_103	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	103	0.62344843
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_104	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	104	0.63793296
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_105	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	105	0.65241754
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_106	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	106	0.6669021
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_107	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	107	0.68138665
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_108	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	108	0.69587123
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_109	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	109	0.71035576
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_110	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	110	0.72484034
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_111	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	111	0.73932487
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_112	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	112	0.75380945
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_113	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	113	0.768294
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_114	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	114	0.78277856
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_115	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	115	0.79726315
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_116	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	116	0.81174767
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_117	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	117	0.82623225
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_118	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	118	0.8407168
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_119	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	119	0.85520136
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_120	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	120	0.8696859
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_121	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	121	0.8841705
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_122	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	122	0.89865506
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_123	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	123	0.9131396
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_124	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	124	0.92762417
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_125	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	125	0.9421087
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_126	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	126	0.9565933
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_127	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	127	0.9710778
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_128	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	128	0.9855624
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_129	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	129	1.000047
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_130	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	130	1.0145315
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_131	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	131	1.029016
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_132	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	132	1.0435007
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_133	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	133	1.0579852
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_134	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	134	1.0724697
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_135	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	135	1.0869542
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_136	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	136	1.1014389
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_137	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	137	1.1159234
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_138	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	138	1.1304079
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_139	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	139	1.1448925
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_141	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	141	1.1738616
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_142	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	142	1.1883461
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_143	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	143	1.2028308
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_144	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	144	1.2173153
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_145	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	145	1.2317998
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_146	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	146	1.2462844
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_147	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	147	1.260769
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_148	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	148	1.2752535
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_149	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	149	1.289738
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_150	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	150	1.3042227
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_151	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	151	1.3187072
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_152	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	152	1.3331918
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_153	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	153	1.3476763
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_154	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	154	1.3621609
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_155	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	155	1.3766454
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_156	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	156	1.39113
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_157	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	157	1.403696
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_158	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	158	1.416262
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_159	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	159	1.428828
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_160	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	160	1.441394
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_161	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	161	1.45396
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_162	Ni2MnGa_sample1_meas_flux_over_field_2	flux_density	Ni2MnGa_sample_1	162	1.466526
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_28	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	28	0.73539776
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_29	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	29	0.7588167
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_30	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	30	0.78223574
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_31	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	31	0.80565476
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_1	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	1	0
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_2	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	2	0.11283355
mechanic_stress_sample1_elastomer1_tens_1_51	tens_1	mechanic_stress	sample1_elastomer1	51	243300
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_3	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	3	0.14992249
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_4	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	4	0.1733415
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_5	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	5	0.1967605
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_6	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	6	0.22017951
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_7	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	7	0.24359852
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_8	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	8	0.26701754
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_9	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	9	0.29043654
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_10	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	10	0.31385556
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_11	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	11	0.33727455
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_12	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	12	0.36069357
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_13	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	13	0.3841126
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_14	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	14	0.4075316
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_15	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	15	0.4309506
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_16	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	16	0.4543696
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_17	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	17	0.47778863
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_18	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	18	0.50120765
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_19	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	19	0.5246266
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_20	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	20	0.54804564
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_21	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	21	0.57146466
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_22	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	22	0.5948837
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_23	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	23	0.6183027
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_24	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	24	0.64172167
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_25	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	25	0.6651407
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_26	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	26	0.6885597
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_27	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	27	0.71197873
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_32	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	32	0.8290738
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_33	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	33	0.85249275
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_34	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	34	0.8759118
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_35	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	35	0.8993308
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_36	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	36	0.9227498
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_37	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	37	0.94616884
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_38	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	38	0.9695878
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_39	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	39	0.9930068
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_40	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	40	1.0164258
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_41	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	41	1.0398449
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_42	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	42	1.0632639
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_43	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	43	1.0866829
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_44	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	44	1.1101019
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_45	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	45	1.1335208
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_46	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	46	1.1569399
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_47	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	47	1.1803589
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_48	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	48	1.2037779
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_49	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	49	1.2271969
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_50	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	50	1.250616
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_51	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	51	1.274035
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_52	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	52	1.297454
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_53	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	53	1.320873
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_54	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	54	1.3442919
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_55	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	55	1.367711
mechanic_stress_sample1_elastomer1_tens_1_52	tens_1	mechanic_stress	sample1_elastomer1	52	247300
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_56	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	56	1.39113
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_57	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	57	1.403696
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_58	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	58	1.416262
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_59	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	59	1.428828
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_60	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	60	1.441394
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_61	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	61	1.45396
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_62	Ni2MnGa_sample1_meas_flux_over_field_1	flux_density	Ni2MnGa_sample_1	62	1.466526
mechanic_stress_sample1_elastomer1_tens_1_1	tens_1	mechanic_stress	sample1_elastomer1	1	0
mechanic_stress_sample1_elastomer1_tens_1_2	tens_1	mechanic_stress	sample1_elastomer1	2	5600
mechanic_stress_sample1_elastomer1_tens_1_3	tens_1	mechanic_stress	sample1_elastomer1	3	12150
mechanic_stress_sample1_elastomer1_tens_1_4	tens_1	mechanic_stress	sample1_elastomer1	4	17200
mechanic_stress_sample1_elastomer1_tens_1_5	tens_1	mechanic_stress	sample1_elastomer1	5	24567
mechanic_stress_sample1_elastomer1_tens_1_6	tens_1	mechanic_stress	sample1_elastomer1	6	30400
mechanic_stress_sample1_elastomer1_tens_1_7	tens_1	mechanic_stress	sample1_elastomer1	7	35900
mechanic_stress_sample1_elastomer1_tens_1_8	tens_1	mechanic_stress	sample1_elastomer1	8	41533
mechanic_stress_sample1_elastomer1_tens_1_9	tens_1	mechanic_stress	sample1_elastomer1	9	47333
mechanic_stress_sample1_elastomer1_tens_1_10	tens_1	mechanic_stress	sample1_elastomer1	10	52900
mechanic_stress_sample1_elastomer1_tens_1_11	tens_1	mechanic_stress	sample1_elastomer1	11	59400
mechanic_stress_sample1_elastomer1_tens_1_12	tens_1	mechanic_stress	sample1_elastomer1	12	64900
mechanic_stress_sample1_elastomer1_tens_1_13	tens_1	mechanic_stress	sample1_elastomer1	13	71000
mechanic_stress_sample1_elastomer1_tens_1_14	tens_1	mechanic_stress	sample1_elastomer1	14	76033
mechanic_stress_sample1_elastomer1_tens_1_15	tens_1	mechanic_stress	sample1_elastomer1	15	82300
mechanic_stress_sample1_elastomer1_tens_1_16	tens_1	mechanic_stress	sample1_elastomer1	16	87300
mechanic_stress_sample1_elastomer1_tens_1_17	tens_1	mechanic_stress	sample1_elastomer1	17	92300
mechanic_stress_sample1_elastomer1_tens_1_18	tens_1	mechanic_stress	sample1_elastomer1	18	97967
mechanic_stress_sample1_elastomer1_tens_1_19	tens_1	mechanic_stress	sample1_elastomer1	19	103300
mechanic_stress_sample1_elastomer1_tens_1_20	tens_1	mechanic_stress	sample1_elastomer1	20	108633
mechanic_stress_sample1_elastomer1_tens_1_21	tens_1	mechanic_stress	sample1_elastomer1	21	113300
mechanic_stress_sample1_elastomer1_tens_1_22	tens_1	mechanic_stress	sample1_elastomer1	22	118300
mechanic_stress_sample1_elastomer1_tens_1_23	tens_1	mechanic_stress	sample1_elastomer1	23	123300
mechanic_stress_sample1_elastomer1_tens_1_24	tens_1	mechanic_stress	sample1_elastomer1	24	129300
mechanic_stress_sample1_elastomer1_tens_1_25	tens_1	mechanic_stress	sample1_elastomer1	25	133300
mechanic_stress_sample1_elastomer1_tens_1_26	tens_1	mechanic_stress	sample1_elastomer1	26	138300
mechanic_stress_sample1_elastomer1_tens_1_27	tens_1	mechanic_stress	sample1_elastomer1	27	143300
mechanic_stress_sample1_elastomer1_tens_1_28	tens_1	mechanic_stress	sample1_elastomer1	28	147300
mechanic_stress_sample1_elastomer1_tens_1_29	tens_1	mechanic_stress	sample1_elastomer1	29	151300
mechanic_stress_sample1_elastomer1_tens_1_30	tens_1	mechanic_stress	sample1_elastomer1	30	156300
mechanic_stress_sample1_elastomer1_tens_1_31	tens_1	mechanic_stress	sample1_elastomer1	31	161300
mechanic_stress_sample1_elastomer1_tens_1_32	tens_1	mechanic_stress	sample1_elastomer1	32	166300
mechanic_stress_sample1_elastomer1_tens_1_33	tens_1	mechanic_stress	sample1_elastomer1	33	170300
mechanic_stress_sample1_elastomer1_tens_1_34	tens_1	mechanic_stress	sample1_elastomer1	34	174300
mechanic_stress_sample1_elastomer1_tens_1_35	tens_1	mechanic_stress	sample1_elastomer1	35	178300
mechanic_stress_sample1_elastomer1_tens_1_36	tens_1	mechanic_stress	sample1_elastomer1	36	183300
mechanic_stress_sample1_elastomer1_tens_1_37	tens_1	mechanic_stress	sample1_elastomer1	37	187300
mechanic_stress_sample1_elastomer1_tens_1_38	tens_1	mechanic_stress	sample1_elastomer1	38	191300
mechanic_stress_sample1_elastomer1_tens_1_39	tens_1	mechanic_stress	sample1_elastomer1	39	196300
mechanic_stress_sample1_elastomer1_tens_1_40	tens_1	mechanic_stress	sample1_elastomer1	40	200300
mechanic_stress_sample1_elastomer1_tens_1_41	tens_1	mechanic_stress	sample1_elastomer1	41	204300
mechanic_stress_sample1_elastomer1_tens_1_42	tens_1	mechanic_stress	sample1_elastomer1	42	208300
mechanic_stress_sample1_elastomer1_tens_1_43	tens_1	mechanic_stress	sample1_elastomer1	43	213300
mechanic_stress_sample1_elastomer1_tens_1_44	tens_1	mechanic_stress	sample1_elastomer1	44	216300
mechanic_stress_sample1_elastomer1_tens_1_45	tens_1	mechanic_stress	sample1_elastomer1	45	220300
mechanic_stress_sample1_elastomer1_tens_1_46	tens_1	mechanic_stress	sample1_elastomer1	46	224300
mechanic_stress_sample1_elastomer1_tens_1_47	tens_1	mechanic_stress	sample1_elastomer1	47	229300
mechanic_stress_sample1_elastomer1_tens_1_48	tens_1	mechanic_stress	sample1_elastomer1	48	231300
mechanic_stress_sample1_elastomer1_tens_1_49	tens_1	mechanic_stress	sample1_elastomer1	49	235300
mechanic_stress_sample1_elastomer1_tens_1_50	tens_1	mechanic_stress	sample1_elastomer1	50	241300
mechanic_stress_sample1_elastomer1_tens_1_53	tens_1	mechanic_stress	sample1_elastomer1	53	250300
mechanic_stress_sample1_elastomer1_tens_1_54	tens_1	mechanic_stress	sample1_elastomer1	54	254300
temperature_wire_actuator_1_123	dummy-meas	temperature	wire_actuator_1	123	75
mechanic_stress_sample1_elastomer1_tens_1_55	tens_1	mechanic_stress	sample1_elastomer1	55	257300
mechanic_stress_sample1_elastomer1_tens_1_56	tens_1	mechanic_stress	sample1_elastomer1	56	262300
mechanic_stress_sample1_elastomer1_tens_1_57	tens_1	mechanic_stress	sample1_elastomer1	57	265300
mechanic_stress_sample1_elastomer1_tens_1_58	tens_1	mechanic_stress	sample1_elastomer1	58	269300
mechanic_stress_sample1_elastomer1_tens_1_59	tens_1	mechanic_stress	sample1_elastomer1	59	273300
mechanic_stress_sample1_elastomer1_tens_1_60	tens_1	mechanic_stress	sample1_elastomer1	60	277300
mechanic_stress_sample1_elastomer1_tens_1_61	tens_1	mechanic_stress	sample1_elastomer1	61	281300
mechanic_stress_sample1_elastomer1_tens_1_62	tens_1	mechanic_stress	sample1_elastomer1	62	283300
mechanic_stress_sample1_elastomer1_tens_1_63	tens_1	mechanic_stress	sample1_elastomer1	63	287300
mechanic_stress_sample1_elastomer1_tens_1_64	tens_1	mechanic_stress	sample1_elastomer1	64	290300
mechanic_stress_sample1_elastomer1_tens_1_65	tens_1	mechanic_stress	sample1_elastomer1	65	294300
mechanic_stress_sample1_elastomer1_tens_1_66	tens_1	mechanic_stress	sample1_elastomer1	66	297300
mechanic_stress_sample1_elastomer1_tens_1_67	tens_1	mechanic_stress	sample1_elastomer1	67	300300
mechanic_stress_sample1_elastomer1_tens_1_68	tens_1	mechanic_stress	sample1_elastomer1	68	303300
mechanic_stress_sample1_elastomer1_tens_1_69	tens_1	mechanic_stress	sample1_elastomer1	69	307300
mechanic_stress_sample1_elastomer1_tens_1_70	tens_1	mechanic_stress	sample1_elastomer1	70	310300
mechanic_stress_sample1_elastomer1_tens_1_71	tens_1	mechanic_stress	sample1_elastomer1	71	313300
mechanic_stress_sample1_elastomer1_tens_1_72	tens_1	mechanic_stress	sample1_elastomer1	72	316300
mechanic_stress_sample1_elastomer1_tens_1_73	tens_1	mechanic_stress	sample1_elastomer1	73	319300
mechanic_stress_sample1_elastomer1_tens_1_74	tens_1	mechanic_stress	sample1_elastomer1	74	323300
mechanic_stress_sample1_elastomer1_tens_1_75	tens_1	mechanic_stress	sample1_elastomer1	75	327300
mechanic_stress_sample1_elastomer1_tens_1_76	tens_1	mechanic_stress	sample1_elastomer1	76	329300
mechanic_stress_sample1_elastomer1_tens_1_77	tens_1	mechanic_stress	sample1_elastomer1	77	333300
mechanic_stress_sample1_elastomer1_tens_1_78	tens_1	mechanic_stress	sample1_elastomer1	78	337300
mechanic_stress_sample1_elastomer1_tens_1_79	tens_1	mechanic_stress	sample1_elastomer1	79	339300
mechanic_stress_sample1_elastomer1_tens_1_80	tens_1	mechanic_stress	sample1_elastomer1	80	342300
mechanic_stress_sample1_elastomer1_tens_1_81	tens_1	mechanic_stress	sample1_elastomer1	81	344300
mechanic_stress_sample1_elastomer1_tens_1_82	tens_1	mechanic_stress	sample1_elastomer1	82	348300
mechanic_stress_sample1_elastomer1_tens_1_83	tens_1	mechanic_stress	sample1_elastomer1	83	351300
mechanic_stress_sample1_elastomer1_tens_1_84	tens_1	mechanic_stress	sample1_elastomer1	84	353300
mechanic_stress_sample1_elastomer1_tens_1_85	tens_1	mechanic_stress	sample1_elastomer1	85	356300
mechanic_stress_sample1_elastomer1_tens_1_86	tens_1	mechanic_stress	sample1_elastomer1	86	359300
mechanic_stress_sample1_elastomer1_tens_1_87	tens_1	mechanic_stress	sample1_elastomer1	87	362300
mechanic_stress_sample1_elastomer1_tens_1_88	tens_1	mechanic_stress	sample1_elastomer1	88	366300
mechanic_stress_sample1_elastomer1_tens_1_89	tens_1	mechanic_stress	sample1_elastomer1	89	369300
mechanic_stress_sample1_elastomer1_tens_1_90	tens_1	mechanic_stress	sample1_elastomer1	90	371300
mechanic_stress_sample1_elastomer1_tens_1_91	tens_1	mechanic_stress	sample1_elastomer1	91	375300
mechanic_stress_sample1_elastomer1_tens_1_92	tens_1	mechanic_stress	sample1_elastomer1	92	378300
mechanic_stress_sample1_elastomer1_tens_1_93	tens_1	mechanic_stress	sample1_elastomer1	93	380300
mechanic_stress_sample1_elastomer1_tens_1_94	tens_1	mechanic_stress	sample1_elastomer1	94	383300
mechanic_stress_sample1_elastomer1_tens_1_95	tens_1	mechanic_stress	sample1_elastomer1	95	385300
mechanic_stress_sample1_elastomer1_tens_1_96	tens_1	mechanic_stress	sample1_elastomer1	96	388300
mechanic_stress_sample1_elastomer1_tens_1_97	tens_1	mechanic_stress	sample1_elastomer1	97	391300
mechanic_stress_sample1_elastomer1_tens_1_98	tens_1	mechanic_stress	sample1_elastomer1	98	394300
mechanic_stress_sample1_elastomer1_tens_1_99	tens_1	mechanic_stress	sample1_elastomer1	99	397300
mechanic_stress_sample1_elastomer1_tens_1_100	tens_1	mechanic_stress	sample1_elastomer1	100	400300
mechanic_stress_sample1_elastomer1_tens_1_101	tens_1	mechanic_stress	sample1_elastomer1	101	402300
mechanic_stress_sample1_elastomer1_tens_1_102	tens_1	mechanic_stress	sample1_elastomer1	102	406300
mechanic_stress_sample1_elastomer1_tens_1_103	tens_1	mechanic_stress	sample1_elastomer1	103	409300
mechanic_stress_sample1_elastomer1_tens_1_104	tens_1	mechanic_stress	sample1_elastomer1	104	410300
mechanic_stress_sample1_elastomer1_tens_1_105	tens_1	mechanic_stress	sample1_elastomer1	105	414300
mechanic_stress_sample1_elastomer1_tens_1_106	tens_1	mechanic_stress	sample1_elastomer1	106	416300
mechanic_stress_sample1_elastomer1_tens_1_107	tens_1	mechanic_stress	sample1_elastomer1	107	419300
mechanic_stress_sample1_elastomer1_tens_1_108	tens_1	mechanic_stress	sample1_elastomer1	108	422300
mechanic_stress_sample1_elastomer1_tens_1_109	tens_1	mechanic_stress	sample1_elastomer1	109	423300
mechanic_stress_sample1_elastomer1_tens_1_110	tens_1	mechanic_stress	sample1_elastomer1	110	428300
mechanic_stress_sample1_elastomer1_tens_1_111	tens_1	mechanic_stress	sample1_elastomer1	111	430300
mechanic_stress_sample1_elastomer1_tens_1_112	tens_1	mechanic_stress	sample1_elastomer1	112	433300
mechanic_stress_sample1_elastomer1_tens_1_113	tens_1	mechanic_stress	sample1_elastomer1	113	435300
mechanic_stress_sample1_elastomer1_tens_1_114	tens_1	mechanic_stress	sample1_elastomer1	114	438300
mechanic_stress_sample1_elastomer1_tens_1_115	tens_1	mechanic_stress	sample1_elastomer1	115	441300
mechanic_stress_sample1_elastomer1_tens_1_116	tens_1	mechanic_stress	sample1_elastomer1	116	443300
mechanic_stress_sample1_elastomer1_tens_1_117	tens_1	mechanic_stress	sample1_elastomer1	117	446300
mechanic_stress_sample1_elastomer1_tens_1_118	tens_1	mechanic_stress	sample1_elastomer1	118	448300
mechanic_stress_sample1_elastomer1_tens_1_119	tens_1	mechanic_stress	sample1_elastomer1	119	452300
mechanic_stress_sample1_elastomer1_tens_1_120	tens_1	mechanic_stress	sample1_elastomer1	120	454300
mechanic_stress_sample1_elastomer1_tens_1_121	tens_1	mechanic_stress	sample1_elastomer1	121	457300
mechanic_stress_sample1_elastomer1_tens_1_122	tens_1	mechanic_stress	sample1_elastomer1	122	460300
mechanic_stress_sample1_elastomer1_tens_1_123	tens_1	mechanic_stress	sample1_elastomer1	123	463300
mechanic_stress_sample1_elastomer1_tens_1_124	tens_1	mechanic_stress	sample1_elastomer1	124	464300
mechanic_stress_sample1_elastomer1_tens_1_125	tens_1	mechanic_stress	sample1_elastomer1	125	468300
mechanic_stress_sample1_elastomer1_tens_1_126	tens_1	mechanic_stress	sample1_elastomer1	126	471300
mechanic_stress_sample1_elastomer1_tens_1_127	tens_1	mechanic_stress	sample1_elastomer1	127	474300
mechanic_stress_sample1_elastomer1_tens_1_128	tens_1	mechanic_stress	sample1_elastomer1	128	475300
mechanic_stress_sample1_elastomer1_tens_1_129	tens_1	mechanic_stress	sample1_elastomer1	129	477300
mechanic_stress_sample1_elastomer1_tens_1_130	tens_1	mechanic_stress	sample1_elastomer1	130	479300
mechanic_stress_sample1_elastomer1_tens_1_131	tens_1	mechanic_stress	sample1_elastomer1	131	483300
mechanic_stress_sample1_elastomer1_tens_1_132	tens_1	mechanic_stress	sample1_elastomer1	132	485300
mechanic_stress_sample1_elastomer1_tens_1_133	tens_1	mechanic_stress	sample1_elastomer1	133	488300
mechanic_stress_sample1_elastomer1_tens_1_134	tens_1	mechanic_stress	sample1_elastomer1	134	491300
mechanic_stress_sample1_elastomer1_tens_1_135	tens_1	mechanic_stress	sample1_elastomer1	135	493300
mechanic_stress_sample1_elastomer1_tens_1_136	tens_1	mechanic_stress	sample1_elastomer1	136	497300
mechanic_stress_sample1_elastomer1_tens_1_137	tens_1	mechanic_stress	sample1_elastomer1	137	499300
mechanic_stress_sample1_elastomer1_tens_1_138	tens_1	mechanic_stress	sample1_elastomer1	138	502300
mechanic_stress_sample1_elastomer1_tens_1_139	tens_1	mechanic_stress	sample1_elastomer1	139	503300
mechanic_stress_sample1_elastomer1_tens_1_140	tens_1	mechanic_stress	sample1_elastomer1	140	506300
mechanic_stress_sample1_elastomer1_tens_1_141	tens_1	mechanic_stress	sample1_elastomer1	141	509300
mechanic_stress_sample1_elastomer1_tens_1_142	tens_1	mechanic_stress	sample1_elastomer1	142	512300
mechanic_stress_sample1_elastomer1_tens_1_143	tens_1	mechanic_stress	sample1_elastomer1	143	513300
mechanic_stress_sample1_elastomer1_tens_1_144	tens_1	mechanic_stress	sample1_elastomer1	144	517300
mechanic_stress_sample1_elastomer1_tens_1_145	tens_1	mechanic_stress	sample1_elastomer1	145	519300
mechanic_stress_sample1_elastomer1_tens_1_146	tens_1	mechanic_stress	sample1_elastomer1	146	522300
mechanic_stress_sample1_elastomer1_tens_1_147	tens_1	mechanic_stress	sample1_elastomer1	147	525300
mechanic_stress_sample1_elastomer1_tens_1_148	tens_1	mechanic_stress	sample1_elastomer1	148	527300
mechanic_stress_sample1_elastomer1_tens_1_149	tens_1	mechanic_stress	sample1_elastomer1	149	530300
mechanic_stress_sample1_elastomer1_tens_1_150	tens_1	mechanic_stress	sample1_elastomer1	150	532300
mechanic_stress_sample1_elastomer1_tens_1_151	tens_1	mechanic_stress	sample1_elastomer1	151	535300
mechanic_stress_sample1_elastomer1_tens_1_152	tens_1	mechanic_stress	sample1_elastomer1	152	538300
mechanic_stress_sample1_elastomer1_tens_1_153	tens_1	mechanic_stress	sample1_elastomer1	153	539300
mechanic_stress_sample1_elastomer1_tens_1_154	tens_1	mechanic_stress	sample1_elastomer1	154	542300
mechanic_stress_sample1_elastomer1_tens_1_155	tens_1	mechanic_stress	sample1_elastomer1	155	545300
mechanic_stress_sample1_elastomer1_tens_1_156	tens_1	mechanic_stress	sample1_elastomer1	156	547300
mechanic_stress_sample1_elastomer1_tens_1_157	tens_1	mechanic_stress	sample1_elastomer1	157	550300
mechanic_stress_sample1_elastomer1_tens_1_158	tens_1	mechanic_stress	sample1_elastomer1	158	553300
mechanic_stress_sample1_elastomer1_tens_1_159	tens_1	mechanic_stress	sample1_elastomer1	159	555300
mechanic_stress_sample1_elastomer1_tens_1_160	tens_1	mechanic_stress	sample1_elastomer1	160	558300
mechanic_stress_sample1_elastomer1_tens_1_161	tens_1	mechanic_stress	sample1_elastomer1	161	560300
mechanic_stress_sample1_elastomer1_tens_1_162	tens_1	mechanic_stress	sample1_elastomer1	162	562300
mechanic_stress_sample1_elastomer1_tens_1_163	tens_1	mechanic_stress	sample1_elastomer1	163	565300
mechanic_stress_sample1_elastomer1_tens_1_164	tens_1	mechanic_stress	sample1_elastomer1	164	568300
mechanic_stress_sample1_elastomer1_tens_1_165	tens_1	mechanic_stress	sample1_elastomer1	165	570300
mechanic_stress_sample1_elastomer1_tens_1_166	tens_1	mechanic_stress	sample1_elastomer1	166	572300
mechanic_stress_sample1_elastomer1_tens_1_167	tens_1	mechanic_stress	sample1_elastomer1	167	576300
mechanic_stress_sample1_elastomer1_tens_1_168	tens_1	mechanic_stress	sample1_elastomer1	168	578300
mechanic_stress_sample1_elastomer1_tens_1_169	tens_1	mechanic_stress	sample1_elastomer1	169	581300
mechanic_stress_sample1_elastomer1_tens_1_170	tens_1	mechanic_stress	sample1_elastomer1	170	584300
mechanic_stress_sample1_elastomer1_tens_1_171	tens_1	mechanic_stress	sample1_elastomer1	171	586300
mechanic_stress_sample1_elastomer1_tens_1_172	tens_1	mechanic_stress	sample1_elastomer1	172	589300
mechanic_stress_sample1_elastomer1_tens_1_173	tens_1	mechanic_stress	sample1_elastomer1	173	591300
mechanic_stress_sample1_elastomer1_tens_1_174	tens_1	mechanic_stress	sample1_elastomer1	174	594300
mechanic_stress_sample1_elastomer1_tens_1_175	tens_1	mechanic_stress	sample1_elastomer1	175	596300
mechanic_stress_sample1_elastomer1_tens_1_176	tens_1	mechanic_stress	sample1_elastomer1	176	598300
mechanic_stress_sample1_elastomer1_tens_1_177	tens_1	mechanic_stress	sample1_elastomer1	177	602300
mechanic_stress_sample1_elastomer1_tens_1_178	tens_1	mechanic_stress	sample1_elastomer1	178	604300
mechanic_stress_sample1_elastomer1_tens_1_179	tens_1	mechanic_stress	sample1_elastomer1	179	607300
mechanic_stress_sample1_elastomer1_tens_1_180	tens_1	mechanic_stress	sample1_elastomer1	180	610300
mechanic_stress_sample1_elastomer1_tens_1_181	tens_1	mechanic_stress	sample1_elastomer1	181	612300
mechanic_stress_sample1_elastomer1_tens_1_182	tens_1	mechanic_stress	sample1_elastomer1	182	613300
mechanic_stress_sample1_elastomer1_tens_1_183	tens_1	mechanic_stress	sample1_elastomer1	183	618300
mechanic_stress_sample1_elastomer1_tens_1_184	tens_1	mechanic_stress	sample1_elastomer1	184	621300
mechanic_stress_sample1_elastomer1_tens_1_185	tens_1	mechanic_stress	sample1_elastomer1	185	622300
mechanic_stress_sample1_elastomer1_tens_1_186	tens_1	mechanic_stress	sample1_elastomer1	186	625300
mechanic_stress_sample1_elastomer1_tens_1_187	tens_1	mechanic_stress	sample1_elastomer1	187	628300
mechanic_stress_sample1_elastomer1_tens_1_188	tens_1	mechanic_stress	sample1_elastomer1	188	630300
mechanic_stress_sample1_elastomer1_tens_1_189	tens_1	mechanic_stress	sample1_elastomer1	189	633300
mechanic_stress_sample1_elastomer1_tens_1_190	tens_1	mechanic_stress	sample1_elastomer1	190	636300
mechanic_stress_sample1_elastomer1_tens_1_191	tens_1	mechanic_stress	sample1_elastomer1	191	638300
mechanic_stress_sample1_elastomer1_tens_1_192	tens_1	mechanic_stress	sample1_elastomer1	192	641300
mechanic_stress_sample1_elastomer1_tens_1_193	tens_1	mechanic_stress	sample1_elastomer1	193	643300
mechanic_stress_sample1_elastomer1_tens_1_194	tens_1	mechanic_stress	sample1_elastomer1	194	646300
mechanic_stress_sample1_elastomer1_tens_1_195	tens_1	mechanic_stress	sample1_elastomer1	195	649300
mechanic_stress_sample1_elastomer1_tens_1_196	tens_1	mechanic_stress	sample1_elastomer1	196	651300
mechanic_stress_sample1_elastomer1_tens_1_197	tens_1	mechanic_stress	sample1_elastomer1	197	654300
mechanic_stress_sample1_elastomer1_tens_1_198	tens_1	mechanic_stress	sample1_elastomer1	198	656300
mechanic_stress_sample1_elastomer1_tens_1_199	tens_1	mechanic_stress	sample1_elastomer1	199	659300
mechanic_stress_sample1_elastomer1_tens_1_200	tens_1	mechanic_stress	sample1_elastomer1	200	663300
mechanic_stress_sample1_elastomer1_tens_1_201	tens_1	mechanic_stress	sample1_elastomer1	201	665300
mechanic_stress_sample1_elastomer1_tens_1_202	tens_1	mechanic_stress	sample1_elastomer1	202	666800
mechanic_stress_sample1_elastomer1_tens_1_203	tens_1	mechanic_stress	sample1_elastomer1	203	668300
mechanic_stress_sample1_elastomer1_tens_1_204	tens_1	mechanic_stress	sample1_elastomer1	204	670800
mechanic_stress_sample1_elastomer1_tens_1_205	tens_1	mechanic_stress	sample1_elastomer1	205	673300
mechanic_stress_sample1_elastomer1_tens_1_206	tens_1	mechanic_stress	sample1_elastomer1	206	676300
mechanic_stress_sample1_elastomer1_tens_1_207	tens_1	mechanic_stress	sample1_elastomer1	207	679300
mechanic_stress_sample1_elastomer1_tens_1_208	tens_1	mechanic_stress	sample1_elastomer1	208	682300
mechanic_stress_sample1_elastomer1_tens_1_209	tens_1	mechanic_stress	sample1_elastomer1	209	685300
mechanic_stress_sample1_elastomer1_tens_1_210	tens_1	mechanic_stress	sample1_elastomer1	210	688300
mechanic_stress_sample1_elastomer1_tens_1_211	tens_1	mechanic_stress	sample1_elastomer1	211	691300
mechanic_stress_sample1_elastomer1_tens_1_212	tens_1	mechanic_stress	sample1_elastomer1	212	693800
mechanic_stress_sample1_elastomer1_tens_1_213	tens_1	mechanic_stress	sample1_elastomer1	213	696300
mechanic_stress_sample1_elastomer1_tens_1_214	tens_1	mechanic_stress	sample1_elastomer1	214	698800
mechanic_stress_sample1_elastomer1_tens_1_215	tens_1	mechanic_stress	sample1_elastomer1	215	701300
mechanic_stress_sample1_elastomer1_tens_1_216	tens_1	mechanic_stress	sample1_elastomer1	216	704300
mechanic_stress_sample1_elastomer1_tens_1_217	tens_1	mechanic_stress	sample1_elastomer1	217	707300
mechanic_stress_sample1_elastomer1_tens_1_218	tens_1	mechanic_stress	sample1_elastomer1	218	709800
mechanic_stress_sample1_elastomer1_tens_1_219	tens_1	mechanic_stress	sample1_elastomer1	219	712300
mechanic_stress_sample1_elastomer1_tens_1_220	tens_1	mechanic_stress	sample1_elastomer1	220	715300
mechanic_stress_sample1_elastomer1_tens_1_221	tens_1	mechanic_stress	sample1_elastomer1	221	718300
mechanic_stress_sample1_elastomer1_tens_1_222	tens_1	mechanic_stress	sample1_elastomer1	222	720800
mechanic_stress_sample1_elastomer1_tens_1_223	tens_1	mechanic_stress	sample1_elastomer1	223	723300
mechanic_stress_sample1_elastomer1_tens_1_224	tens_1	mechanic_stress	sample1_elastomer1	224	726300
mechanic_stress_sample1_elastomer1_tens_1_225	tens_1	mechanic_stress	sample1_elastomer1	225	729300
mechanic_stress_sample1_elastomer1_tens_1_226	tens_1	mechanic_stress	sample1_elastomer1	226	732300
mechanic_stress_sample1_elastomer1_tens_1_227	tens_1	mechanic_stress	sample1_elastomer1	227	735300
mechanic_stress_sample1_elastomer1_tens_1_228	tens_1	mechanic_stress	sample1_elastomer1	228	737800
mechanic_stress_sample1_elastomer1_tens_1_229	tens_1	mechanic_stress	sample1_elastomer1	229	740300
mechanic_stress_sample1_elastomer1_tens_1_230	tens_1	mechanic_stress	sample1_elastomer1	230	742800
mechanic_stress_sample1_elastomer1_tens_1_231	tens_1	mechanic_stress	sample1_elastomer1	231	745300
mechanic_stress_sample1_elastomer1_tens_1_232	tens_1	mechanic_stress	sample1_elastomer1	232	747800
mechanic_stress_sample1_elastomer1_tens_1_233	tens_1	mechanic_stress	sample1_elastomer1	233	750300
mechanic_stress_sample1_elastomer1_tens_1_234	tens_1	mechanic_stress	sample1_elastomer1	234	753300
mechanic_stress_sample1_elastomer1_tens_1_235	tens_1	mechanic_stress	sample1_elastomer1	235	756300
mechanic_stress_sample1_elastomer1_tens_1_236	tens_1	mechanic_stress	sample1_elastomer1	236	759800
mechanic_stress_sample1_elastomer1_tens_1_237	tens_1	mechanic_stress	sample1_elastomer1	237	763300
mechanic_stress_sample1_elastomer1_tens_1_238	tens_1	mechanic_stress	sample1_elastomer1	238	765800
mechanic_stress_sample1_elastomer1_tens_1_239	tens_1	mechanic_stress	sample1_elastomer1	239	768300
mechanic_stress_sample1_elastomer1_tens_1_240	tens_1	mechanic_stress	sample1_elastomer1	240	771300
mechanic_stress_sample1_elastomer1_tens_1_241	tens_1	mechanic_stress	sample1_elastomer1	241	774300
mechanic_stress_sample1_elastomer1_tens_1_242	tens_1	mechanic_stress	sample1_elastomer1	242	776300
mechanic_stress_sample1_elastomer1_tens_1_243	tens_1	mechanic_stress	sample1_elastomer1	243	778300
mechanic_stress_sample1_elastomer1_tens_1_244	tens_1	mechanic_stress	sample1_elastomer1	244	782300
mechanic_stress_sample1_elastomer1_tens_1_245	tens_1	mechanic_stress	sample1_elastomer1	245	786300
mechanic_stress_sample1_elastomer1_tens_1_246	tens_1	mechanic_stress	sample1_elastomer1	246	789300
mechanic_stress_sample1_elastomer1_tens_1_247	tens_1	mechanic_stress	sample1_elastomer1	247	792300
mechanic_stress_sample1_elastomer1_tens_1_248	tens_1	mechanic_stress	sample1_elastomer1	248	795300
mechanic_stress_sample1_elastomer1_tens_1_249	tens_1	mechanic_stress	sample1_elastomer1	249	798300
mechanic_stress_sample1_elastomer1_tens_1_250	tens_1	mechanic_stress	sample1_elastomer1	250	800800
mechanic_stress_sample1_elastomer1_tens_1_251	tens_1	mechanic_stress	sample1_elastomer1	251	803300
mechanic_stress_sample1_elastomer1_tens_1_252	tens_1	mechanic_stress	sample1_elastomer1	252	806800
mechanic_stress_sample1_elastomer1_tens_1_253	tens_1	mechanic_stress	sample1_elastomer1	253	810300
mechanic_stress_sample1_elastomer1_tens_1_254	tens_1	mechanic_stress	sample1_elastomer1	254	812800
mechanic_stress_sample1_elastomer1_tens_1_255	tens_1	mechanic_stress	sample1_elastomer1	255	815300
mechanic_stress_sample1_elastomer1_tens_1_256	tens_1	mechanic_stress	sample1_elastomer1	256	818300
mechanic_stress_sample1_elastomer1_tens_1_257	tens_1	mechanic_stress	sample1_elastomer1	257	821300
mechanic_stress_sample1_elastomer1_tens_1_258	tens_1	mechanic_stress	sample1_elastomer1	258	824300
mechanic_stress_sample1_elastomer1_tens_1_259	tens_1	mechanic_stress	sample1_elastomer1	259	827300
mechanic_stress_sample1_elastomer1_tens_1_260	tens_1	mechanic_stress	sample1_elastomer1	260	830300
mechanic_stress_sample1_elastomer1_tens_1_261	tens_1	mechanic_stress	sample1_elastomer1	261	833300
mechanic_stress_sample1_elastomer1_tens_1_262	tens_1	mechanic_stress	sample1_elastomer1	262	836300
mechanic_stress_sample1_elastomer1_tens_1_263	tens_1	mechanic_stress	sample1_elastomer1	263	839300
mechanic_stress_sample1_elastomer1_tens_1_264	tens_1	mechanic_stress	sample1_elastomer1	264	842300
mechanic_stress_sample1_elastomer1_tens_1_265	tens_1	mechanic_stress	sample1_elastomer1	265	845300
mechanic_stress_sample1_elastomer1_tens_1_266	tens_1	mechanic_stress	sample1_elastomer1	266	848300
mechanic_stress_sample1_elastomer1_tens_1_267	tens_1	mechanic_stress	sample1_elastomer1	267	851300
mechanic_stress_sample1_elastomer1_tens_1_268	tens_1	mechanic_stress	sample1_elastomer1	268	854300
mechanic_stress_sample1_elastomer1_tens_1_269	tens_1	mechanic_stress	sample1_elastomer1	269	857300
mechanic_stress_sample1_elastomer1_tens_1_270	tens_1	mechanic_stress	sample1_elastomer1	270	860800
mechanic_stress_sample1_elastomer1_tens_1_271	tens_1	mechanic_stress	sample1_elastomer1	271	864300
mechanic_stress_sample1_elastomer1_tens_1_272	tens_1	mechanic_stress	sample1_elastomer1	272	867800
mechanic_stress_sample1_elastomer1_tens_1_273	tens_1	mechanic_stress	sample1_elastomer1	273	871300
mechanic_stress_sample1_elastomer1_tens_1_274	tens_1	mechanic_stress	sample1_elastomer1	274	873800
mechanic_stress_sample1_elastomer1_tens_1_275	tens_1	mechanic_stress	sample1_elastomer1	275	876300
mechanic_stress_sample1_elastomer1_tens_1_276	tens_1	mechanic_stress	sample1_elastomer1	276	878800
mechanic_stress_sample1_elastomer1_tens_1_277	tens_1	mechanic_stress	sample1_elastomer1	277	881300
mechanic_stress_sample1_elastomer1_tens_1_278	tens_1	mechanic_stress	sample1_elastomer1	278	885300
mechanic_stress_sample1_elastomer1_tens_1_279	tens_1	mechanic_stress	sample1_elastomer1	279	889300
mechanic_stress_sample1_elastomer1_tens_1_280	tens_1	mechanic_stress	sample1_elastomer1	280	891800
mechanic_stress_sample1_elastomer1_tens_1_281	tens_1	mechanic_stress	sample1_elastomer1	281	894300
mechanic_stress_sample1_elastomer1_tens_1_282	tens_1	mechanic_stress	sample1_elastomer1	282	897800
mechanic_stress_sample1_elastomer1_tens_1_283	tens_1	mechanic_stress	sample1_elastomer1	283	901300
mechanic_stress_sample1_elastomer1_tens_1_284	tens_1	mechanic_stress	sample1_elastomer1	284	904300
mechanic_stress_sample1_elastomer1_tens_1_285	tens_1	mechanic_stress	sample1_elastomer1	285	907300
mechanic_stress_sample1_elastomer1_tens_1_286	tens_1	mechanic_stress	sample1_elastomer1	286	910300
mechanic_stress_sample1_elastomer1_tens_1_287	tens_1	mechanic_stress	sample1_elastomer1	287	913300
mechanic_stress_sample1_elastomer1_tens_1_288	tens_1	mechanic_stress	sample1_elastomer1	288	916800
mechanic_stress_sample1_elastomer1_tens_1_289	tens_1	mechanic_stress	sample1_elastomer1	289	920300
mechanic_stress_sample1_elastomer1_tens_1_290	tens_1	mechanic_stress	sample1_elastomer1	290	923800
mechanic_stress_sample1_elastomer1_tens_1_291	tens_1	mechanic_stress	sample1_elastomer1	291	927300
mechanic_stress_sample1_elastomer1_tens_1_292	tens_1	mechanic_stress	sample1_elastomer1	292	930300
mechanic_stress_sample1_elastomer1_tens_1_293	tens_1	mechanic_stress	sample1_elastomer1	293	933300
mechanic_stress_sample1_elastomer1_tens_1_294	tens_1	mechanic_stress	sample1_elastomer1	294	936300
mechanic_stress_sample1_elastomer1_tens_1_295	tens_1	mechanic_stress	sample1_elastomer1	295	939300
mechanic_stress_sample1_elastomer1_tens_1_296	tens_1	mechanic_stress	sample1_elastomer1	296	942300
mechanic_stress_sample1_elastomer1_tens_1_297	tens_1	mechanic_stress	sample1_elastomer1	297	945300
mechanic_stress_sample1_elastomer1_tens_1_298	tens_1	mechanic_stress	sample1_elastomer1	298	948800
mechanic_stress_sample1_elastomer1_tens_1_299	tens_1	mechanic_stress	sample1_elastomer1	299	952300
mechanic_stress_sample1_elastomer1_tens_1_300	tens_1	mechanic_stress	sample1_elastomer1	300	955800
mechanic_stress_sample1_elastomer1_tens_1_301	tens_1	mechanic_stress	sample1_elastomer1	301	959300
mechanic_stress_sample1_elastomer1_tens_1_302	tens_1	mechanic_stress	sample1_elastomer1	302	962300
mechanic_stress_sample1_elastomer1_tens_1_303	tens_1	mechanic_stress	sample1_elastomer1	303	965300
mechanic_stress_sample1_elastomer1_tens_1_304	tens_1	mechanic_stress	sample1_elastomer1	304	968800
mechanic_stress_sample1_elastomer1_tens_1_305	tens_1	mechanic_stress	sample1_elastomer1	305	972300
mechanic_stress_sample1_elastomer1_tens_1_306	tens_1	mechanic_stress	sample1_elastomer1	306	975300
mechanic_stress_sample1_elastomer1_tens_1_307	tens_1	mechanic_stress	sample1_elastomer1	307	978300
mechanic_stress_sample1_elastomer1_tens_1_308	tens_1	mechanic_stress	sample1_elastomer1	308	983300
mechanic_stress_sample1_elastomer1_tens_1_309	tens_1	mechanic_stress	sample1_elastomer1	309	988300
mechanic_stress_sample1_elastomer1_tens_1_310	tens_1	mechanic_stress	sample1_elastomer1	310	988300
mechanic_stress_sample1_elastomer1_tens_1_311	tens_1	mechanic_stress	sample1_elastomer1	311	988300
mechanic_stress_sample1_elastomer1_tens_1_312	tens_1	mechanic_stress	sample1_elastomer1	312	993300
mechanic_stress_sample1_elastomer1_tens_1_313	tens_1	mechanic_stress	sample1_elastomer1	313	998300
mechanic_stress_sample1_elastomer1_tens_1_314	tens_1	mechanic_stress	sample1_elastomer1	314	1003300
mechanic_stress_sample1_elastomer1_tens_1_315	tens_1	mechanic_stress	sample1_elastomer1	315	1008300
mechanic_stress_sample1_elastomer1_tens_1_316	tens_1	mechanic_stress	sample1_elastomer1	316	1008300
mechanic_stress_sample1_elastomer1_tens_1_317	tens_1	mechanic_stress	sample1_elastomer1	317	1008300
mechanic_stress_sample1_elastomer1_tens_1_318	tens_1	mechanic_stress	sample1_elastomer1	318	1013300
mechanic_stress_sample1_elastomer1_tens_1_319	tens_1	mechanic_stress	sample1_elastomer1	319	1018300
mechanic_stress_sample1_elastomer1_tens_1_320	tens_1	mechanic_stress	sample1_elastomer1	320	1023300
mechanic_stress_sample1_elastomer1_tens_1_321	tens_1	mechanic_stress	sample1_elastomer1	321	1028300
mechanic_stress_sample1_elastomer1_tens_1_322	tens_1	mechanic_stress	sample1_elastomer1	322	1028300
mechanic_stress_sample1_elastomer1_tens_1_323	tens_1	mechanic_stress	sample1_elastomer1	323	1028300
mechanic_stress_sample1_elastomer1_tens_1_324	tens_1	mechanic_stress	sample1_elastomer1	324	1033300
mechanic_stress_sample1_elastomer1_tens_1_325	tens_1	mechanic_stress	sample1_elastomer1	325	1038300
mechanic_stress_sample1_elastomer1_tens_1_326	tens_1	mechanic_stress	sample1_elastomer1	326	1043300
mechanic_stress_sample1_elastomer1_tens_1_327	tens_1	mechanic_stress	sample1_elastomer1	327	1048300
mechanic_stress_sample1_elastomer1_tens_1_328	tens_1	mechanic_stress	sample1_elastomer1	328	1048300
mechanic_stress_sample1_elastomer1_tens_1_329	tens_1	mechanic_stress	sample1_elastomer1	329	1048300
mechanic_stress_sample1_elastomer1_tens_1_330	tens_1	mechanic_stress	sample1_elastomer1	330	1053300
mechanic_stress_sample1_elastomer1_tens_1_331	tens_1	mechanic_stress	sample1_elastomer1	331	1058300
mechanic_stress_sample1_elastomer1_tens_1_332	tens_1	mechanic_stress	sample1_elastomer1	332	1063300
mechanic_stress_sample1_elastomer1_tens_1_333	tens_1	mechanic_stress	sample1_elastomer1	333	1068300
mechanic_stress_sample1_elastomer1_tens_1_334	tens_1	mechanic_stress	sample1_elastomer1	334	1073300
mechanic_stress_sample1_elastomer1_tens_1_335	tens_1	mechanic_stress	sample1_elastomer1	335	1078300
mechanic_stress_sample1_elastomer1_tens_1_336	tens_1	mechanic_stress	sample1_elastomer1	336	1078300
mechanic_stress_sample1_elastomer1_tens_1_337	tens_1	mechanic_stress	sample1_elastomer1	337	1078300
mechanic_stress_sample1_elastomer1_tens_1_338	tens_1	mechanic_stress	sample1_elastomer1	338	1083300
mechanic_stress_sample1_elastomer1_tens_1_339	tens_1	mechanic_stress	sample1_elastomer1	339	1088300
mechanic_stress_sample1_elastomer1_tens_1_340	tens_1	mechanic_stress	sample1_elastomer1	340	1088300
mechanic_stress_sample1_elastomer1_tens_1_341	tens_1	mechanic_stress	sample1_elastomer1	341	1088300
mechanic_stress_sample1_elastomer1_tens_1_342	tens_1	mechanic_stress	sample1_elastomer1	342	1093300
mechanic_stress_sample1_elastomer1_tens_1_343	tens_1	mechanic_stress	sample1_elastomer1	343	1098300
mechanic_stress_sample1_elastomer1_tens_1_344	tens_1	mechanic_stress	sample1_elastomer1	344	1103300
mechanic_stress_sample1_elastomer1_tens_1_345	tens_1	mechanic_stress	sample1_elastomer1	345	1108300
mechanic_stress_sample1_elastomer1_tens_1_346	tens_1	mechanic_stress	sample1_elastomer1	346	1113300
mechanic_stress_sample1_elastomer1_tens_1_347	tens_1	mechanic_stress	sample1_elastomer1	347	1118300
mechanic_stress_sample1_elastomer1_tens_1_348	tens_1	mechanic_stress	sample1_elastomer1	348	1118300
mechanic_stress_sample1_elastomer1_tens_1_349	tens_1	mechanic_stress	sample1_elastomer1	349	1118300
mechanic_stress_sample1_elastomer1_tens_1_350	tens_1	mechanic_stress	sample1_elastomer1	350	1123300
mechanic_stress_sample1_elastomer1_tens_1_351	tens_1	mechanic_stress	sample1_elastomer1	351	1128300
mechanic_stress_sample1_elastomer1_tens_1_352	tens_1	mechanic_stress	sample1_elastomer1	352	1128300
mechanic_stress_sample1_elastomer1_tens_1_353	tens_1	mechanic_stress	sample1_elastomer1	353	1128300
mechanic_stress_sample1_elastomer1_tens_1_354	tens_1	mechanic_stress	sample1_elastomer1	354	1133300
mechanic_stress_sample1_elastomer1_tens_1_355	tens_1	mechanic_stress	sample1_elastomer1	355	1138300
mechanic_stress_sample1_elastomer1_tens_1_356	tens_1	mechanic_stress	sample1_elastomer1	356	1143300
mechanic_stress_sample1_elastomer1_tens_1_357	tens_1	mechanic_stress	sample1_elastomer1	357	1148300
mechanic_stress_sample1_elastomer1_tens_1_358	tens_1	mechanic_stress	sample1_elastomer1	358	1153300
mechanic_stress_sample1_elastomer1_tens_1_359	tens_1	mechanic_stress	sample1_elastomer1	359	1158300
mechanic_stress_sample1_elastomer1_tens_1_360	tens_1	mechanic_stress	sample1_elastomer1	360	1158300
mechanic_stress_sample1_elastomer1_tens_1_361	tens_1	mechanic_stress	sample1_elastomer1	361	1158300
mechanic_stress_sample1_elastomer1_tens_1_362	tens_1	mechanic_stress	sample1_elastomer1	362	1163300
mechanic_stress_sample1_elastomer1_tens_1_363	tens_1	mechanic_stress	sample1_elastomer1	363	1168300
mechanic_stress_sample1_elastomer1_tens_1_364	tens_1	mechanic_stress	sample1_elastomer1	364	1173300
mechanic_stress_sample1_elastomer1_tens_1_365	tens_1	mechanic_stress	sample1_elastomer1	365	1178300
mechanic_stress_sample1_elastomer1_tens_1_366	tens_1	mechanic_stress	sample1_elastomer1	366	1178300
mechanic_stress_sample1_elastomer1_tens_1_367	tens_1	mechanic_stress	sample1_elastomer1	367	1178300
mechanic_stress_sample1_elastomer1_tens_1_368	tens_1	mechanic_stress	sample1_elastomer1	368	1183300
mechanic_stress_sample1_elastomer1_tens_1_369	tens_1	mechanic_stress	sample1_elastomer1	369	1188300
mechanic_stress_sample1_elastomer1_tens_1_370	tens_1	mechanic_stress	sample1_elastomer1	370	1193300
mechanic_stress_sample1_elastomer1_tens_1_371	tens_1	mechanic_stress	sample1_elastomer1	371	1198300
mechanic_stress_sample1_elastomer1_tens_1_372	tens_1	mechanic_stress	sample1_elastomer1	372	1198300
mechanic_stress_sample1_elastomer1_tens_1_373	tens_1	mechanic_stress	sample1_elastomer1	373	1198300
mechanic_stress_sample1_elastomer1_tens_1_374	tens_1	mechanic_stress	sample1_elastomer1	374	1203300
mechanic_stress_sample1_elastomer1_tens_1_375	tens_1	mechanic_stress	sample1_elastomer1	375	1208300
mechanic_stress_sample1_elastomer1_tens_1_376	tens_1	mechanic_stress	sample1_elastomer1	376	1213300
mechanic_stress_sample1_elastomer1_tens_1_377	tens_1	mechanic_stress	sample1_elastomer1	377	1218300
mechanic_stress_sample1_elastomer1_tens_1_378	tens_1	mechanic_stress	sample1_elastomer1	378	1218300
mechanic_stress_sample1_elastomer1_tens_1_379	tens_1	mechanic_stress	sample1_elastomer1	379	1218300
mechanic_stress_sample1_elastomer1_tens_1_380	tens_1	mechanic_stress	sample1_elastomer1	380	1223300
mechanic_stress_sample1_elastomer1_tens_1_381	tens_1	mechanic_stress	sample1_elastomer1	381	1228300
mechanic_stress_sample1_elastomer1_tens_1_382	tens_1	mechanic_stress	sample1_elastomer1	382	1233300
mechanic_stress_sample1_elastomer1_tens_1_383	tens_1	mechanic_stress	sample1_elastomer1	383	1238300
mechanic_stress_sample1_elastomer1_tens_1_384	tens_1	mechanic_stress	sample1_elastomer1	384	1243300
mechanic_stress_sample1_elastomer1_tens_1_385	tens_1	mechanic_stress	sample1_elastomer1	385	1248300
mechanic_stress_sample1_elastomer1_tens_1_386	tens_1	mechanic_stress	sample1_elastomer1	386	1248300
mechanic_stress_sample1_elastomer1_tens_1_387	tens_1	mechanic_stress	sample1_elastomer1	387	1248300
mechanic_stress_sample1_elastomer1_tens_1_388	tens_1	mechanic_stress	sample1_elastomer1	388	1253300
mechanic_stress_sample1_elastomer1_tens_1_389	tens_1	mechanic_stress	sample1_elastomer1	389	1258300
mechanic_stress_sample1_elastomer1_tens_1_390	tens_1	mechanic_stress	sample1_elastomer1	390	1263300
mechanic_stress_sample1_elastomer1_tens_1_391	tens_1	mechanic_stress	sample1_elastomer1	391	1268300
mechanic_stress_sample1_elastomer1_tens_1_392	tens_1	mechanic_stress	sample1_elastomer1	392	1268300
mechanic_stress_sample1_elastomer1_tens_1_393	tens_1	mechanic_stress	sample1_elastomer1	393	1268300
mechanic_stress_sample1_elastomer1_tens_1_394	tens_1	mechanic_stress	sample1_elastomer1	394	1273300
mechanic_stress_sample1_elastomer1_tens_1_395	tens_1	mechanic_stress	sample1_elastomer1	395	1278300
mechanic_stress_sample1_elastomer1_tens_1_396	tens_1	mechanic_stress	sample1_elastomer1	396	1283300
mechanic_stress_sample1_elastomer1_tens_1_397	tens_1	mechanic_stress	sample1_elastomer1	397	1288300
mechanic_stress_sample1_elastomer1_tens_1_398	tens_1	mechanic_stress	sample1_elastomer1	398	1293300
mechanic_stress_sample1_elastomer1_tens_1_399	tens_1	mechanic_stress	sample1_elastomer1	399	1298300
mechanic_stress_sample1_elastomer1_tens_1_400	tens_1	mechanic_stress	sample1_elastomer1	400	1303300
mechanic_stress_sample1_elastomer1_tens_1_401	tens_1	mechanic_stress	sample1_elastomer1	401	1308300
strain_sample1_elastomer1_tens_1_24	tens_1	strain	sample1_elastomer1	24	0.115
strain_sample1_elastomer1_tens_1_1	tens_1	strain	sample1_elastomer1	1	0
strain_sample1_elastomer1_tens_1_2	tens_1	strain	sample1_elastomer1	2	0.005
strain_sample1_elastomer1_tens_1_3	tens_1	strain	sample1_elastomer1	3	0.01
strain_sample1_elastomer1_tens_1_4	tens_1	strain	sample1_elastomer1	4	0.015
strain_sample1_elastomer1_tens_1_5	tens_1	strain	sample1_elastomer1	5	0.02
strain_sample1_elastomer1_tens_1_6	tens_1	strain	sample1_elastomer1	6	0.025
strain_sample1_elastomer1_tens_1_7	tens_1	strain	sample1_elastomer1	7	0.03
strain_sample1_elastomer1_tens_1_8	tens_1	strain	sample1_elastomer1	8	0.035
strain_sample1_elastomer1_tens_1_9	tens_1	strain	sample1_elastomer1	9	0.04
strain_sample1_elastomer1_tens_1_10	tens_1	strain	sample1_elastomer1	10	0.045
strain_sample1_elastomer1_tens_1_11	tens_1	strain	sample1_elastomer1	11	0.05
strain_sample1_elastomer1_tens_1_12	tens_1	strain	sample1_elastomer1	12	0.055
strain_sample1_elastomer1_tens_1_13	tens_1	strain	sample1_elastomer1	13	0.06
strain_sample1_elastomer1_tens_1_14	tens_1	strain	sample1_elastomer1	14	0.065
strain_sample1_elastomer1_tens_1_15	tens_1	strain	sample1_elastomer1	15	0.07
strain_sample1_elastomer1_tens_1_16	tens_1	strain	sample1_elastomer1	16	0.075
strain_sample1_elastomer1_tens_1_17	tens_1	strain	sample1_elastomer1	17	0.08
strain_sample1_elastomer1_tens_1_18	tens_1	strain	sample1_elastomer1	18	0.085
strain_sample1_elastomer1_tens_1_19	tens_1	strain	sample1_elastomer1	19	0.09
strain_sample1_elastomer1_tens_1_20	tens_1	strain	sample1_elastomer1	20	0.095
strain_sample1_elastomer1_tens_1_21	tens_1	strain	sample1_elastomer1	21	0.1
strain_sample1_elastomer1_tens_1_22	tens_1	strain	sample1_elastomer1	22	0.105
strain_sample1_elastomer1_tens_1_23	tens_1	strain	sample1_elastomer1	23	0.11
strain_sample1_elastomer1_tens_1_25	tens_1	strain	sample1_elastomer1	25	0.12
strain_sample1_elastomer1_tens_1_26	tens_1	strain	sample1_elastomer1	26	0.125
strain_sample1_elastomer1_tens_1_27	tens_1	strain	sample1_elastomer1	27	0.13
strain_sample1_elastomer1_tens_1_28	tens_1	strain	sample1_elastomer1	28	0.135
strain_sample1_elastomer1_tens_1_29	tens_1	strain	sample1_elastomer1	29	0.14
strain_sample1_elastomer1_tens_1_30	tens_1	strain	sample1_elastomer1	30	0.145
strain_sample1_elastomer1_tens_1_31	tens_1	strain	sample1_elastomer1	31	0.15
strain_sample1_elastomer1_tens_1_32	tens_1	strain	sample1_elastomer1	32	0.155
strain_sample1_elastomer1_tens_1_33	tens_1	strain	sample1_elastomer1	33	0.16
strain_sample1_elastomer1_tens_1_34	tens_1	strain	sample1_elastomer1	34	0.165
strain_sample1_elastomer1_tens_1_35	tens_1	strain	sample1_elastomer1	35	0.17
strain_sample1_elastomer1_tens_1_36	tens_1	strain	sample1_elastomer1	36	0.175
strain_sample1_elastomer1_tens_1_37	tens_1	strain	sample1_elastomer1	37	0.18
strain_sample1_elastomer1_tens_1_38	tens_1	strain	sample1_elastomer1	38	0.185
strain_sample1_elastomer1_tens_1_39	tens_1	strain	sample1_elastomer1	39	0.19
strain_sample1_elastomer1_tens_1_40	tens_1	strain	sample1_elastomer1	40	0.195
strain_sample1_elastomer1_tens_1_41	tens_1	strain	sample1_elastomer1	41	0.2
strain_sample1_elastomer1_tens_1_42	tens_1	strain	sample1_elastomer1	42	0.205
strain_sample1_elastomer1_tens_1_43	tens_1	strain	sample1_elastomer1	43	0.21
strain_sample1_elastomer1_tens_1_44	tens_1	strain	sample1_elastomer1	44	0.215
strain_sample1_elastomer1_tens_1_45	tens_1	strain	sample1_elastomer1	45	0.22
strain_sample1_elastomer1_tens_1_46	tens_1	strain	sample1_elastomer1	46	0.225
strain_sample1_elastomer1_tens_1_47	tens_1	strain	sample1_elastomer1	47	0.23
strain_sample1_elastomer1_tens_1_48	tens_1	strain	sample1_elastomer1	48	0.235
strain_sample1_elastomer1_tens_1_49	tens_1	strain	sample1_elastomer1	49	0.24
strain_sample1_elastomer1_tens_1_50	tens_1	strain	sample1_elastomer1	50	0.245
strain_sample1_elastomer1_tens_1_51	tens_1	strain	sample1_elastomer1	51	0.25
strain_sample1_elastomer1_tens_1_52	tens_1	strain	sample1_elastomer1	52	0.255
strain_sample1_elastomer1_tens_1_53	tens_1	strain	sample1_elastomer1	53	0.26
strain_sample1_elastomer1_tens_1_54	tens_1	strain	sample1_elastomer1	54	0.265
strain_sample1_elastomer1_tens_1_55	tens_1	strain	sample1_elastomer1	55	0.27
strain_sample1_elastomer1_tens_1_56	tens_1	strain	sample1_elastomer1	56	0.275
strain_sample1_elastomer1_tens_1_57	tens_1	strain	sample1_elastomer1	57	0.28
strain_sample1_elastomer1_tens_1_58	tens_1	strain	sample1_elastomer1	58	0.285
strain_sample1_elastomer1_tens_1_59	tens_1	strain	sample1_elastomer1	59	0.29
strain_sample1_elastomer1_tens_1_60	tens_1	strain	sample1_elastomer1	60	0.295
strain_sample1_elastomer1_tens_1_61	tens_1	strain	sample1_elastomer1	61	0.3
strain_sample1_elastomer1_tens_1_62	tens_1	strain	sample1_elastomer1	62	0.305
strain_sample1_elastomer1_tens_1_63	tens_1	strain	sample1_elastomer1	63	0.31
strain_sample1_elastomer1_tens_1_64	tens_1	strain	sample1_elastomer1	64	0.315
strain_sample1_elastomer1_tens_1_65	tens_1	strain	sample1_elastomer1	65	0.32
strain_sample1_elastomer1_tens_1_66	tens_1	strain	sample1_elastomer1	66	0.325
strain_sample1_elastomer1_tens_1_67	tens_1	strain	sample1_elastomer1	67	0.33
strain_sample1_elastomer1_tens_1_68	tens_1	strain	sample1_elastomer1	68	0.335
strain_sample1_elastomer1_tens_1_69	tens_1	strain	sample1_elastomer1	69	0.34
strain_sample1_elastomer1_tens_1_70	tens_1	strain	sample1_elastomer1	70	0.345
strain_sample1_elastomer1_tens_1_71	tens_1	strain	sample1_elastomer1	71	0.35
strain_sample1_elastomer1_tens_1_72	tens_1	strain	sample1_elastomer1	72	0.355
strain_sample1_elastomer1_tens_1_73	tens_1	strain	sample1_elastomer1	73	0.36
strain_sample1_elastomer1_tens_1_74	tens_1	strain	sample1_elastomer1	74	0.365
strain_sample1_elastomer1_tens_1_75	tens_1	strain	sample1_elastomer1	75	0.37
strain_sample1_elastomer1_tens_1_76	tens_1	strain	sample1_elastomer1	76	0.375
strain_sample1_elastomer1_tens_1_77	tens_1	strain	sample1_elastomer1	77	0.38
strain_sample1_elastomer1_tens_1_78	tens_1	strain	sample1_elastomer1	78	0.385
strain_sample1_elastomer1_tens_1_79	tens_1	strain	sample1_elastomer1	79	0.39
strain_sample1_elastomer1_tens_1_80	tens_1	strain	sample1_elastomer1	80	0.395
strain_sample1_elastomer1_tens_1_81	tens_1	strain	sample1_elastomer1	81	0.4
strain_sample1_elastomer1_tens_1_82	tens_1	strain	sample1_elastomer1	82	0.405
strain_sample1_elastomer1_tens_1_83	tens_1	strain	sample1_elastomer1	83	0.41
strain_sample1_elastomer1_tens_1_84	tens_1	strain	sample1_elastomer1	84	0.415
strain_sample1_elastomer1_tens_1_85	tens_1	strain	sample1_elastomer1	85	0.42
strain_sample1_elastomer1_tens_1_86	tens_1	strain	sample1_elastomer1	86	0.425
strain_sample1_elastomer1_tens_1_87	tens_1	strain	sample1_elastomer1	87	0.43
strain_sample1_elastomer1_tens_1_88	tens_1	strain	sample1_elastomer1	88	0.435
strain_sample1_elastomer1_tens_1_89	tens_1	strain	sample1_elastomer1	89	0.44
strain_sample1_elastomer1_tens_1_90	tens_1	strain	sample1_elastomer1	90	0.445
strain_sample1_elastomer1_tens_1_91	tens_1	strain	sample1_elastomer1	91	0.45
strain_sample1_elastomer1_tens_1_92	tens_1	strain	sample1_elastomer1	92	0.455
strain_sample1_elastomer1_tens_1_93	tens_1	strain	sample1_elastomer1	93	0.46
strain_sample1_elastomer1_tens_1_94	tens_1	strain	sample1_elastomer1	94	0.465
strain_sample1_elastomer1_tens_1_95	tens_1	strain	sample1_elastomer1	95	0.47
strain_sample1_elastomer1_tens_1_96	tens_1	strain	sample1_elastomer1	96	0.475
strain_sample1_elastomer1_tens_1_97	tens_1	strain	sample1_elastomer1	97	0.48
strain_sample1_elastomer1_tens_1_98	tens_1	strain	sample1_elastomer1	98	0.485
strain_sample1_elastomer1_tens_1_99	tens_1	strain	sample1_elastomer1	99	0.49
strain_sample1_elastomer1_tens_1_100	tens_1	strain	sample1_elastomer1	100	0.495
strain_sample1_elastomer1_tens_1_101	tens_1	strain	sample1_elastomer1	101	0.5
strain_sample1_elastomer1_tens_1_102	tens_1	strain	sample1_elastomer1	102	0.505
strain_sample1_elastomer1_tens_1_103	tens_1	strain	sample1_elastomer1	103	0.51
strain_sample1_elastomer1_tens_1_104	tens_1	strain	sample1_elastomer1	104	0.515
strain_sample1_elastomer1_tens_1_105	tens_1	strain	sample1_elastomer1	105	0.52
strain_sample1_elastomer1_tens_1_106	tens_1	strain	sample1_elastomer1	106	0.525
strain_sample1_elastomer1_tens_1_107	tens_1	strain	sample1_elastomer1	107	0.53
strain_sample1_elastomer1_tens_1_108	tens_1	strain	sample1_elastomer1	108	0.535
strain_sample1_elastomer1_tens_1_109	tens_1	strain	sample1_elastomer1	109	0.54
strain_sample1_elastomer1_tens_1_110	tens_1	strain	sample1_elastomer1	110	0.545
strain_sample1_elastomer1_tens_1_111	tens_1	strain	sample1_elastomer1	111	0.55
strain_sample1_elastomer1_tens_1_112	tens_1	strain	sample1_elastomer1	112	0.555
strain_sample1_elastomer1_tens_1_113	tens_1	strain	sample1_elastomer1	113	0.56
strain_sample1_elastomer1_tens_1_114	tens_1	strain	sample1_elastomer1	114	0.565
strain_sample1_elastomer1_tens_1_115	tens_1	strain	sample1_elastomer1	115	0.57
strain_sample1_elastomer1_tens_1_116	tens_1	strain	sample1_elastomer1	116	0.575
strain_sample1_elastomer1_tens_1_117	tens_1	strain	sample1_elastomer1	117	0.58
strain_sample1_elastomer1_tens_1_118	tens_1	strain	sample1_elastomer1	118	0.585
strain_sample1_elastomer1_tens_1_119	tens_1	strain	sample1_elastomer1	119	0.59
strain_sample1_elastomer1_tens_1_120	tens_1	strain	sample1_elastomer1	120	0.595
strain_sample1_elastomer1_tens_1_121	tens_1	strain	sample1_elastomer1	121	0.6
strain_sample1_elastomer1_tens_1_122	tens_1	strain	sample1_elastomer1	122	0.605
strain_sample1_elastomer1_tens_1_123	tens_1	strain	sample1_elastomer1	123	0.61
strain_sample1_elastomer1_tens_1_124	tens_1	strain	sample1_elastomer1	124	0.615
strain_sample1_elastomer1_tens_1_125	tens_1	strain	sample1_elastomer1	125	0.62
strain_sample1_elastomer1_tens_1_126	tens_1	strain	sample1_elastomer1	126	0.625
strain_sample1_elastomer1_tens_1_127	tens_1	strain	sample1_elastomer1	127	0.63
strain_sample1_elastomer1_tens_1_128	tens_1	strain	sample1_elastomer1	128	0.635
strain_sample1_elastomer1_tens_1_129	tens_1	strain	sample1_elastomer1	129	0.64
strain_sample1_elastomer1_tens_1_130	tens_1	strain	sample1_elastomer1	130	0.645
strain_sample1_elastomer1_tens_1_131	tens_1	strain	sample1_elastomer1	131	0.65
strain_sample1_elastomer1_tens_1_132	tens_1	strain	sample1_elastomer1	132	0.655
strain_sample1_elastomer1_tens_1_133	tens_1	strain	sample1_elastomer1	133	0.66
strain_sample1_elastomer1_tens_1_134	tens_1	strain	sample1_elastomer1	134	0.665
strain_sample1_elastomer1_tens_1_135	tens_1	strain	sample1_elastomer1	135	0.67
strain_sample1_elastomer1_tens_1_136	tens_1	strain	sample1_elastomer1	136	0.675
strain_sample1_elastomer1_tens_1_137	tens_1	strain	sample1_elastomer1	137	0.68
strain_sample1_elastomer1_tens_1_138	tens_1	strain	sample1_elastomer1	138	0.685
strain_sample1_elastomer1_tens_1_139	tens_1	strain	sample1_elastomer1	139	0.69
strain_sample1_elastomer1_tens_1_140	tens_1	strain	sample1_elastomer1	140	0.695
strain_sample1_elastomer1_tens_1_141	tens_1	strain	sample1_elastomer1	141	0.7
strain_sample1_elastomer1_tens_1_142	tens_1	strain	sample1_elastomer1	142	0.705
strain_sample1_elastomer1_tens_1_143	tens_1	strain	sample1_elastomer1	143	0.71
strain_sample1_elastomer1_tens_1_144	tens_1	strain	sample1_elastomer1	144	0.715
strain_sample1_elastomer1_tens_1_145	tens_1	strain	sample1_elastomer1	145	0.72
strain_sample1_elastomer1_tens_1_146	tens_1	strain	sample1_elastomer1	146	0.725
strain_sample1_elastomer1_tens_1_147	tens_1	strain	sample1_elastomer1	147	0.73
strain_sample1_elastomer1_tens_1_148	tens_1	strain	sample1_elastomer1	148	0.735
strain_sample1_elastomer1_tens_1_149	tens_1	strain	sample1_elastomer1	149	0.74
strain_sample1_elastomer1_tens_1_150	tens_1	strain	sample1_elastomer1	150	0.745
strain_sample1_elastomer1_tens_1_151	tens_1	strain	sample1_elastomer1	151	0.75
strain_sample1_elastomer1_tens_1_152	tens_1	strain	sample1_elastomer1	152	0.755
strain_sample1_elastomer1_tens_1_153	tens_1	strain	sample1_elastomer1	153	0.76
strain_sample1_elastomer1_tens_1_154	tens_1	strain	sample1_elastomer1	154	0.765
strain_sample1_elastomer1_tens_1_155	tens_1	strain	sample1_elastomer1	155	0.77
strain_sample1_elastomer1_tens_1_156	tens_1	strain	sample1_elastomer1	156	0.775
strain_sample1_elastomer1_tens_1_157	tens_1	strain	sample1_elastomer1	157	0.78
strain_sample1_elastomer1_tens_1_158	tens_1	strain	sample1_elastomer1	158	0.785
strain_sample1_elastomer1_tens_1_159	tens_1	strain	sample1_elastomer1	159	0.79
strain_sample1_elastomer1_tens_1_160	tens_1	strain	sample1_elastomer1	160	0.795
strain_sample1_elastomer1_tens_1_161	tens_1	strain	sample1_elastomer1	161	0.8
strain_sample1_elastomer1_tens_1_162	tens_1	strain	sample1_elastomer1	162	0.805
strain_sample1_elastomer1_tens_1_163	tens_1	strain	sample1_elastomer1	163	0.81
strain_sample1_elastomer1_tens_1_164	tens_1	strain	sample1_elastomer1	164	0.815
strain_sample1_elastomer1_tens_1_165	tens_1	strain	sample1_elastomer1	165	0.82
strain_sample1_elastomer1_tens_1_166	tens_1	strain	sample1_elastomer1	166	0.825
strain_sample1_elastomer1_tens_1_167	tens_1	strain	sample1_elastomer1	167	0.83
strain_sample1_elastomer1_tens_1_168	tens_1	strain	sample1_elastomer1	168	0.835
strain_sample1_elastomer1_tens_1_169	tens_1	strain	sample1_elastomer1	169	0.84
strain_sample1_elastomer1_tens_1_170	tens_1	strain	sample1_elastomer1	170	0.845
strain_sample1_elastomer1_tens_1_171	tens_1	strain	sample1_elastomer1	171	0.85
strain_sample1_elastomer1_tens_1_172	tens_1	strain	sample1_elastomer1	172	0.855
strain_sample1_elastomer1_tens_1_173	tens_1	strain	sample1_elastomer1	173	0.86
strain_sample1_elastomer1_tens_1_174	tens_1	strain	sample1_elastomer1	174	0.865
strain_sample1_elastomer1_tens_1_175	tens_1	strain	sample1_elastomer1	175	0.87
strain_sample1_elastomer1_tens_1_176	tens_1	strain	sample1_elastomer1	176	0.875
strain_sample1_elastomer1_tens_1_177	tens_1	strain	sample1_elastomer1	177	0.88
strain_sample1_elastomer1_tens_1_178	tens_1	strain	sample1_elastomer1	178	0.885
strain_sample1_elastomer1_tens_1_179	tens_1	strain	sample1_elastomer1	179	0.89
strain_sample1_elastomer1_tens_1_180	tens_1	strain	sample1_elastomer1	180	0.895
strain_sample1_elastomer1_tens_1_181	tens_1	strain	sample1_elastomer1	181	0.9
strain_sample1_elastomer1_tens_1_182	tens_1	strain	sample1_elastomer1	182	0.905
strain_sample1_elastomer1_tens_1_183	tens_1	strain	sample1_elastomer1	183	0.91
strain_sample1_elastomer1_tens_1_184	tens_1	strain	sample1_elastomer1	184	0.915
strain_sample1_elastomer1_tens_1_185	tens_1	strain	sample1_elastomer1	185	0.92
strain_sample1_elastomer1_tens_1_186	tens_1	strain	sample1_elastomer1	186	0.925
strain_sample1_elastomer1_tens_1_187	tens_1	strain	sample1_elastomer1	187	0.93
strain_sample1_elastomer1_tens_1_188	tens_1	strain	sample1_elastomer1	188	0.935
strain_sample1_elastomer1_tens_1_189	tens_1	strain	sample1_elastomer1	189	0.94
strain_sample1_elastomer1_tens_1_190	tens_1	strain	sample1_elastomer1	190	0.945
strain_sample1_elastomer1_tens_1_191	tens_1	strain	sample1_elastomer1	191	0.95
strain_sample1_elastomer1_tens_1_192	tens_1	strain	sample1_elastomer1	192	0.955
strain_sample1_elastomer1_tens_1_193	tens_1	strain	sample1_elastomer1	193	0.96
strain_sample1_elastomer1_tens_1_194	tens_1	strain	sample1_elastomer1	194	0.965
strain_sample1_elastomer1_tens_1_195	tens_1	strain	sample1_elastomer1	195	0.97
strain_sample1_elastomer1_tens_1_196	tens_1	strain	sample1_elastomer1	196	0.975
strain_sample1_elastomer1_tens_1_197	tens_1	strain	sample1_elastomer1	197	0.98
strain_sample1_elastomer1_tens_1_198	tens_1	strain	sample1_elastomer1	198	0.985
strain_sample1_elastomer1_tens_1_199	tens_1	strain	sample1_elastomer1	199	0.99
strain_sample1_elastomer1_tens_1_200	tens_1	strain	sample1_elastomer1	200	0.995
strain_sample1_elastomer1_tens_1_201	tens_1	strain	sample1_elastomer1	201	1
strain_sample1_elastomer1_tens_1_202	tens_1	strain	sample1_elastomer1	202	1.005
strain_sample1_elastomer1_tens_1_203	tens_1	strain	sample1_elastomer1	203	1.01
strain_sample1_elastomer1_tens_1_204	tens_1	strain	sample1_elastomer1	204	1.015
strain_sample1_elastomer1_tens_1_205	tens_1	strain	sample1_elastomer1	205	1.02
strain_sample1_elastomer1_tens_1_206	tens_1	strain	sample1_elastomer1	206	1.025
strain_sample1_elastomer1_tens_1_207	tens_1	strain	sample1_elastomer1	207	1.03
strain_sample1_elastomer1_tens_1_208	tens_1	strain	sample1_elastomer1	208	1.035
strain_sample1_elastomer1_tens_1_209	tens_1	strain	sample1_elastomer1	209	1.04
strain_sample1_elastomer1_tens_1_210	tens_1	strain	sample1_elastomer1	210	1.045
strain_sample1_elastomer1_tens_1_211	tens_1	strain	sample1_elastomer1	211	1.05
strain_sample1_elastomer1_tens_1_212	tens_1	strain	sample1_elastomer1	212	1.055
strain_sample1_elastomer1_tens_1_213	tens_1	strain	sample1_elastomer1	213	1.06
strain_sample1_elastomer1_tens_1_214	tens_1	strain	sample1_elastomer1	214	1.065
strain_sample1_elastomer1_tens_1_215	tens_1	strain	sample1_elastomer1	215	1.07
strain_sample1_elastomer1_tens_1_216	tens_1	strain	sample1_elastomer1	216	1.075
strain_sample1_elastomer1_tens_1_217	tens_1	strain	sample1_elastomer1	217	1.08
strain_sample1_elastomer1_tens_1_218	tens_1	strain	sample1_elastomer1	218	1.085
strain_sample1_elastomer1_tens_1_219	tens_1	strain	sample1_elastomer1	219	1.09
strain_sample1_elastomer1_tens_1_220	tens_1	strain	sample1_elastomer1	220	1.095
strain_sample1_elastomer1_tens_1_221	tens_1	strain	sample1_elastomer1	221	1.1
strain_sample1_elastomer1_tens_1_222	tens_1	strain	sample1_elastomer1	222	1.105
strain_sample1_elastomer1_tens_1_223	tens_1	strain	sample1_elastomer1	223	1.11
strain_sample1_elastomer1_tens_1_224	tens_1	strain	sample1_elastomer1	224	1.115
strain_sample1_elastomer1_tens_1_225	tens_1	strain	sample1_elastomer1	225	1.12
strain_sample1_elastomer1_tens_1_226	tens_1	strain	sample1_elastomer1	226	1.125
strain_sample1_elastomer1_tens_1_227	tens_1	strain	sample1_elastomer1	227	1.13
strain_sample1_elastomer1_tens_1_228	tens_1	strain	sample1_elastomer1	228	1.135
strain_sample1_elastomer1_tens_1_229	tens_1	strain	sample1_elastomer1	229	1.14
strain_sample1_elastomer1_tens_1_230	tens_1	strain	sample1_elastomer1	230	1.145
strain_sample1_elastomer1_tens_1_231	tens_1	strain	sample1_elastomer1	231	1.15
strain_sample1_elastomer1_tens_1_232	tens_1	strain	sample1_elastomer1	232	1.155
strain_sample1_elastomer1_tens_1_233	tens_1	strain	sample1_elastomer1	233	1.16
strain_sample1_elastomer1_tens_1_234	tens_1	strain	sample1_elastomer1	234	1.165
strain_sample1_elastomer1_tens_1_235	tens_1	strain	sample1_elastomer1	235	1.17
strain_sample1_elastomer1_tens_1_236	tens_1	strain	sample1_elastomer1	236	1.175
strain_sample1_elastomer1_tens_1_237	tens_1	strain	sample1_elastomer1	237	1.18
strain_sample1_elastomer1_tens_1_238	tens_1	strain	sample1_elastomer1	238	1.185
strain_sample1_elastomer1_tens_1_239	tens_1	strain	sample1_elastomer1	239	1.19
strain_sample1_elastomer1_tens_1_240	tens_1	strain	sample1_elastomer1	240	1.195
strain_sample1_elastomer1_tens_1_241	tens_1	strain	sample1_elastomer1	241	1.2
strain_sample1_elastomer1_tens_1_242	tens_1	strain	sample1_elastomer1	242	1.205
strain_sample1_elastomer1_tens_1_243	tens_1	strain	sample1_elastomer1	243	1.21
strain_sample1_elastomer1_tens_1_244	tens_1	strain	sample1_elastomer1	244	1.215
strain_sample1_elastomer1_tens_1_245	tens_1	strain	sample1_elastomer1	245	1.22
strain_sample1_elastomer1_tens_1_246	tens_1	strain	sample1_elastomer1	246	1.225
strain_sample1_elastomer1_tens_1_247	tens_1	strain	sample1_elastomer1	247	1.23
strain_sample1_elastomer1_tens_1_248	tens_1	strain	sample1_elastomer1	248	1.235
strain_sample1_elastomer1_tens_1_249	tens_1	strain	sample1_elastomer1	249	1.24
strain_sample1_elastomer1_tens_1_250	tens_1	strain	sample1_elastomer1	250	1.245
strain_sample1_elastomer1_tens_1_251	tens_1	strain	sample1_elastomer1	251	1.25
strain_sample1_elastomer1_tens_1_252	tens_1	strain	sample1_elastomer1	252	1.255
strain_sample1_elastomer1_tens_1_253	tens_1	strain	sample1_elastomer1	253	1.26
strain_sample1_elastomer1_tens_1_254	tens_1	strain	sample1_elastomer1	254	1.265
strain_sample1_elastomer1_tens_1_255	tens_1	strain	sample1_elastomer1	255	1.27
strain_sample1_elastomer1_tens_1_256	tens_1	strain	sample1_elastomer1	256	1.275
strain_sample1_elastomer1_tens_1_257	tens_1	strain	sample1_elastomer1	257	1.28
strain_sample1_elastomer1_tens_1_258	tens_1	strain	sample1_elastomer1	258	1.285
strain_sample1_elastomer1_tens_1_259	tens_1	strain	sample1_elastomer1	259	1.29
strain_sample1_elastomer1_tens_1_260	tens_1	strain	sample1_elastomer1	260	1.295
strain_sample1_elastomer1_tens_1_261	tens_1	strain	sample1_elastomer1	261	1.3
strain_sample1_elastomer1_tens_1_262	tens_1	strain	sample1_elastomer1	262	1.305
strain_sample1_elastomer1_tens_1_263	tens_1	strain	sample1_elastomer1	263	1.31
strain_sample1_elastomer1_tens_1_264	tens_1	strain	sample1_elastomer1	264	1.315
strain_sample1_elastomer1_tens_1_265	tens_1	strain	sample1_elastomer1	265	1.32
strain_sample1_elastomer1_tens_1_266	tens_1	strain	sample1_elastomer1	266	1.325
strain_sample1_elastomer1_tens_1_267	tens_1	strain	sample1_elastomer1	267	1.33
strain_sample1_elastomer1_tens_1_268	tens_1	strain	sample1_elastomer1	268	1.335
strain_sample1_elastomer1_tens_1_269	tens_1	strain	sample1_elastomer1	269	1.34
strain_sample1_elastomer1_tens_1_270	tens_1	strain	sample1_elastomer1	270	1.345
strain_sample1_elastomer1_tens_1_271	tens_1	strain	sample1_elastomer1	271	1.35
strain_sample1_elastomer1_tens_1_272	tens_1	strain	sample1_elastomer1	272	1.355
strain_sample1_elastomer1_tens_1_273	tens_1	strain	sample1_elastomer1	273	1.36
strain_sample1_elastomer1_tens_1_274	tens_1	strain	sample1_elastomer1	274	1.365
strain_sample1_elastomer1_tens_1_275	tens_1	strain	sample1_elastomer1	275	1.37
strain_sample1_elastomer1_tens_1_276	tens_1	strain	sample1_elastomer1	276	1.375
strain_sample1_elastomer1_tens_1_277	tens_1	strain	sample1_elastomer1	277	1.38
strain_sample1_elastomer1_tens_1_278	tens_1	strain	sample1_elastomer1	278	1.385
strain_sample1_elastomer1_tens_1_279	tens_1	strain	sample1_elastomer1	279	1.39
strain_sample1_elastomer1_tens_1_280	tens_1	strain	sample1_elastomer1	280	1.395
strain_sample1_elastomer1_tens_1_281	tens_1	strain	sample1_elastomer1	281	1.4
strain_sample1_elastomer1_tens_1_282	tens_1	strain	sample1_elastomer1	282	1.405
strain_sample1_elastomer1_tens_1_283	tens_1	strain	sample1_elastomer1	283	1.41
strain_sample1_elastomer1_tens_1_284	tens_1	strain	sample1_elastomer1	284	1.415
strain_sample1_elastomer1_tens_1_285	tens_1	strain	sample1_elastomer1	285	1.42
strain_sample1_elastomer1_tens_1_286	tens_1	strain	sample1_elastomer1	286	1.425
strain_sample1_elastomer1_tens_1_287	tens_1	strain	sample1_elastomer1	287	1.43
strain_sample1_elastomer1_tens_1_288	tens_1	strain	sample1_elastomer1	288	1.435
strain_sample1_elastomer1_tens_1_289	tens_1	strain	sample1_elastomer1	289	1.44
strain_sample1_elastomer1_tens_1_290	tens_1	strain	sample1_elastomer1	290	1.445
strain_sample1_elastomer1_tens_1_291	tens_1	strain	sample1_elastomer1	291	1.45
strain_sample1_elastomer1_tens_1_292	tens_1	strain	sample1_elastomer1	292	1.455
strain_sample1_elastomer1_tens_1_293	tens_1	strain	sample1_elastomer1	293	1.46
strain_sample1_elastomer1_tens_1_294	tens_1	strain	sample1_elastomer1	294	1.465
strain_sample1_elastomer1_tens_1_295	tens_1	strain	sample1_elastomer1	295	1.47
strain_sample1_elastomer1_tens_1_296	tens_1	strain	sample1_elastomer1	296	1.475
strain_sample1_elastomer1_tens_1_297	tens_1	strain	sample1_elastomer1	297	1.48
strain_sample1_elastomer1_tens_1_298	tens_1	strain	sample1_elastomer1	298	1.485
strain_sample1_elastomer1_tens_1_299	tens_1	strain	sample1_elastomer1	299	1.49
strain_sample1_elastomer1_tens_1_300	tens_1	strain	sample1_elastomer1	300	1.495
strain_sample1_elastomer1_tens_1_301	tens_1	strain	sample1_elastomer1	301	1.5
strain_sample1_elastomer1_tens_1_302	tens_1	strain	sample1_elastomer1	302	1.505
strain_sample1_elastomer1_tens_1_303	tens_1	strain	sample1_elastomer1	303	1.51
strain_sample1_elastomer1_tens_1_304	tens_1	strain	sample1_elastomer1	304	1.515
strain_sample1_elastomer1_tens_1_305	tens_1	strain	sample1_elastomer1	305	1.52
strain_sample1_elastomer1_tens_1_306	tens_1	strain	sample1_elastomer1	306	1.525
strain_sample1_elastomer1_tens_1_307	tens_1	strain	sample1_elastomer1	307	1.53
strain_sample1_elastomer1_tens_1_308	tens_1	strain	sample1_elastomer1	308	1.535
strain_sample1_elastomer1_tens_1_309	tens_1	strain	sample1_elastomer1	309	1.54
strain_sample1_elastomer1_tens_1_310	tens_1	strain	sample1_elastomer1	310	1.545
strain_sample1_elastomer1_tens_1_311	tens_1	strain	sample1_elastomer1	311	1.55
strain_sample1_elastomer1_tens_1_312	tens_1	strain	sample1_elastomer1	312	1.555
strain_sample1_elastomer1_tens_1_313	tens_1	strain	sample1_elastomer1	313	1.56
strain_sample1_elastomer1_tens_1_314	tens_1	strain	sample1_elastomer1	314	1.565
strain_sample1_elastomer1_tens_1_315	tens_1	strain	sample1_elastomer1	315	1.57
strain_sample1_elastomer1_tens_1_316	tens_1	strain	sample1_elastomer1	316	1.575
strain_sample1_elastomer1_tens_1_317	tens_1	strain	sample1_elastomer1	317	1.58
strain_sample1_elastomer1_tens_1_318	tens_1	strain	sample1_elastomer1	318	1.585
strain_sample1_elastomer1_tens_1_319	tens_1	strain	sample1_elastomer1	319	1.59
strain_sample1_elastomer1_tens_1_320	tens_1	strain	sample1_elastomer1	320	1.595
strain_sample1_elastomer1_tens_1_321	tens_1	strain	sample1_elastomer1	321	1.6
strain_sample1_elastomer1_tens_1_322	tens_1	strain	sample1_elastomer1	322	1.605
strain_sample1_elastomer1_tens_1_323	tens_1	strain	sample1_elastomer1	323	1.61
strain_sample1_elastomer1_tens_1_324	tens_1	strain	sample1_elastomer1	324	1.615
strain_sample1_elastomer1_tens_1_325	tens_1	strain	sample1_elastomer1	325	1.62
strain_sample1_elastomer1_tens_1_326	tens_1	strain	sample1_elastomer1	326	1.625
strain_sample1_elastomer1_tens_1_327	tens_1	strain	sample1_elastomer1	327	1.63
strain_sample1_elastomer1_tens_1_328	tens_1	strain	sample1_elastomer1	328	1.635
strain_sample1_elastomer1_tens_1_329	tens_1	strain	sample1_elastomer1	329	1.64
strain_sample1_elastomer1_tens_1_330	tens_1	strain	sample1_elastomer1	330	1.645
strain_sample1_elastomer1_tens_1_331	tens_1	strain	sample1_elastomer1	331	1.65
strain_sample1_elastomer1_tens_1_332	tens_1	strain	sample1_elastomer1	332	1.655
strain_sample1_elastomer1_tens_1_333	tens_1	strain	sample1_elastomer1	333	1.66
strain_sample1_elastomer1_tens_1_334	tens_1	strain	sample1_elastomer1	334	1.665
strain_sample1_elastomer1_tens_1_335	tens_1	strain	sample1_elastomer1	335	1.67
strain_sample1_elastomer1_tens_1_336	tens_1	strain	sample1_elastomer1	336	1.675
strain_sample1_elastomer1_tens_1_337	tens_1	strain	sample1_elastomer1	337	1.68
strain_sample1_elastomer1_tens_1_338	tens_1	strain	sample1_elastomer1	338	1.685
strain_sample1_elastomer1_tens_1_339	tens_1	strain	sample1_elastomer1	339	1.69
strain_sample1_elastomer1_tens_1_340	tens_1	strain	sample1_elastomer1	340	1.695
strain_sample1_elastomer1_tens_1_341	tens_1	strain	sample1_elastomer1	341	1.7
strain_sample1_elastomer1_tens_1_342	tens_1	strain	sample1_elastomer1	342	1.705
strain_sample1_elastomer1_tens_1_343	tens_1	strain	sample1_elastomer1	343	1.71
strain_sample1_elastomer1_tens_1_344	tens_1	strain	sample1_elastomer1	344	1.715
strain_sample1_elastomer1_tens_1_345	tens_1	strain	sample1_elastomer1	345	1.72
strain_sample1_elastomer1_tens_1_346	tens_1	strain	sample1_elastomer1	346	1.725
strain_sample1_elastomer1_tens_1_347	tens_1	strain	sample1_elastomer1	347	1.73
strain_sample1_elastomer1_tens_1_348	tens_1	strain	sample1_elastomer1	348	1.735
strain_sample1_elastomer1_tens_1_349	tens_1	strain	sample1_elastomer1	349	1.74
strain_sample1_elastomer1_tens_1_350	tens_1	strain	sample1_elastomer1	350	1.745
strain_sample1_elastomer1_tens_1_351	tens_1	strain	sample1_elastomer1	351	1.75
strain_sample1_elastomer1_tens_1_352	tens_1	strain	sample1_elastomer1	352	1.755
strain_sample1_elastomer1_tens_1_353	tens_1	strain	sample1_elastomer1	353	1.76
strain_sample1_elastomer1_tens_1_354	tens_1	strain	sample1_elastomer1	354	1.765
strain_sample1_elastomer1_tens_1_355	tens_1	strain	sample1_elastomer1	355	1.77
strain_sample1_elastomer1_tens_1_356	tens_1	strain	sample1_elastomer1	356	1.775
strain_sample1_elastomer1_tens_1_357	tens_1	strain	sample1_elastomer1	357	1.78
strain_sample1_elastomer1_tens_1_358	tens_1	strain	sample1_elastomer1	358	1.785
strain_sample1_elastomer1_tens_1_359	tens_1	strain	sample1_elastomer1	359	1.79
strain_sample1_elastomer1_tens_1_360	tens_1	strain	sample1_elastomer1	360	1.795
strain_sample1_elastomer1_tens_1_361	tens_1	strain	sample1_elastomer1	361	1.8
strain_sample1_elastomer1_tens_1_362	tens_1	strain	sample1_elastomer1	362	1.805
strain_sample1_elastomer1_tens_1_363	tens_1	strain	sample1_elastomer1	363	1.81
strain_sample1_elastomer1_tens_1_364	tens_1	strain	sample1_elastomer1	364	1.815
strain_sample1_elastomer1_tens_1_365	tens_1	strain	sample1_elastomer1	365	1.82
strain_sample1_elastomer1_tens_1_366	tens_1	strain	sample1_elastomer1	366	1.825
strain_sample1_elastomer1_tens_1_367	tens_1	strain	sample1_elastomer1	367	1.83
strain_sample1_elastomer1_tens_1_368	tens_1	strain	sample1_elastomer1	368	1.835
strain_sample1_elastomer1_tens_1_369	tens_1	strain	sample1_elastomer1	369	1.84
strain_sample1_elastomer1_tens_1_370	tens_1	strain	sample1_elastomer1	370	1.845
strain_sample1_elastomer1_tens_1_371	tens_1	strain	sample1_elastomer1	371	1.85
strain_sample1_elastomer1_tens_1_372	tens_1	strain	sample1_elastomer1	372	1.855
strain_sample1_elastomer1_tens_1_373	tens_1	strain	sample1_elastomer1	373	1.86
strain_sample1_elastomer1_tens_1_374	tens_1	strain	sample1_elastomer1	374	1.865
strain_sample1_elastomer1_tens_1_375	tens_1	strain	sample1_elastomer1	375	1.87
strain_sample1_elastomer1_tens_1_376	tens_1	strain	sample1_elastomer1	376	1.875
strain_sample1_elastomer1_tens_1_377	tens_1	strain	sample1_elastomer1	377	1.88
strain_sample1_elastomer1_tens_1_378	tens_1	strain	sample1_elastomer1	378	1.885
strain_sample1_elastomer1_tens_1_379	tens_1	strain	sample1_elastomer1	379	1.89
strain_sample1_elastomer1_tens_1_380	tens_1	strain	sample1_elastomer1	380	1.895
strain_sample1_elastomer1_tens_1_381	tens_1	strain	sample1_elastomer1	381	1.9
strain_sample1_elastomer1_tens_1_382	tens_1	strain	sample1_elastomer1	382	1.905
strain_sample1_elastomer1_tens_1_383	tens_1	strain	sample1_elastomer1	383	1.91
strain_sample1_elastomer1_tens_1_384	tens_1	strain	sample1_elastomer1	384	1.915
strain_sample1_elastomer1_tens_1_385	tens_1	strain	sample1_elastomer1	385	1.92
strain_sample1_elastomer1_tens_1_386	tens_1	strain	sample1_elastomer1	386	1.925
strain_sample1_elastomer1_tens_1_387	tens_1	strain	sample1_elastomer1	387	1.93
strain_sample1_elastomer1_tens_1_388	tens_1	strain	sample1_elastomer1	388	1.935
strain_sample1_elastomer1_tens_1_389	tens_1	strain	sample1_elastomer1	389	1.94
strain_sample1_elastomer1_tens_1_390	tens_1	strain	sample1_elastomer1	390	1.945
strain_sample1_elastomer1_tens_1_391	tens_1	strain	sample1_elastomer1	391	1.95
strain_sample1_elastomer1_tens_1_392	tens_1	strain	sample1_elastomer1	392	1.955
strain_sample1_elastomer1_tens_1_393	tens_1	strain	sample1_elastomer1	393	1.96
strain_sample1_elastomer1_tens_1_394	tens_1	strain	sample1_elastomer1	394	1.965
strain_sample1_elastomer1_tens_1_395	tens_1	strain	sample1_elastomer1	395	1.97
strain_sample1_elastomer1_tens_1_396	tens_1	strain	sample1_elastomer1	396	1.975
strain_sample1_elastomer1_tens_1_397	tens_1	strain	sample1_elastomer1	397	1.98
strain_sample1_elastomer1_tens_1_398	tens_1	strain	sample1_elastomer1	398	1.985
strain_sample1_elastomer1_tens_1_399	tens_1	strain	sample1_elastomer1	399	1.99
strain_sample1_elastomer1_tens_1_400	tens_1	strain	sample1_elastomer1	400	1.995
strain_sample1_elastomer1_tens_1_401	tens_1	strain	sample1_elastomer1	401	2
max_block_stress_M420	concr_calc_max_blocking_stress_PC	maximal_blocking_stress	M420	0	1907689.37
dielectric_strength_el_1	elas1_measured_at	dielectric_strength	elastomer_1	0	99000000
relative_permittivity_el_1	elas1_measured_at	relative_permittivity	elastomer_1	0	3.2
initial_strain_elastomer_1	elas1_base_value	initial_strain	sample1_elastomer1	0	0.05
twinning_stress_stick_Ni2MnGa_sample_1	dummy-meas	twinning_stress	stick_Ni2MnGa_sample_1	0	300000
M420_1_1	dummy-meas	A00	M420	0	1.54e-11
M420_1_2	dummy-meas	A00	M420	0	1.87e-11
M420_1_3	dummy-meas	A00	M420	0	-6.5e-12
M420_1_4	dummy-meas	A00	M420	0	-5.7e-12
temperature_wire_actuator_1_122	dummy-meas	temperature	wire_actuator_1	122	73.5
temperature_wire_actuator_1_121	dummy-meas	temperature	wire_actuator_1	121	72
temperature_wire_actuator_1_120	dummy-meas	temperature	wire_actuator_1	120	71
temperature_wire_actuator_1_119	dummy-meas	temperature	wire_actuator_1	119	70
temperature_wire_actuator_1_118	dummy-meas	temperature	wire_actuator_1	118	65
temperature_wire_actuator_1_117	dummy-meas	temperature	wire_actuator_1	117	60
temperature_wire_actuator_1_116	dummy-meas	temperature	wire_actuator_1	116	53
temperature_wire_actuator_1_115	dummy-meas	temperature	wire_actuator_1	115	52
temperature_wire_actuator_1_111	dummy-meas	temperature	wire_actuator_1	111	45
temperature_wire_actuator_1_114	dummy-meas	temperature	wire_actuator_1	114	51
temperature_wire_actuator_1_113	dummy-meas	temperature	wire_actuator_1	113	50
temperature_wire_actuator_1_112	dummy-meas	temperature	wire_actuator_1	112	48
temperature_wire_actuator_1_131	dummy-meas	temperature	wire_actuator_1	131	45
temperature_wire_actuator_1_140	dummy-meas	temperature	wire_actuator_1	140	71
temperature_wire_actuator_1_141	dummy-meas	temperature	wire_actuator_1	141	72
temperature_wire_actuator_1_142	dummy-meas	temperature	wire_actuator_1	142	73.5
temperature_wire_actuator_1_143	dummy-meas	temperature	wire_actuator_1	143	75
temperature_wire_actuator_1_132	dummy-meas	temperature	wire_actuator_1	132	48
temperature_wire_actuator_1_133	dummy-meas	temperature	wire_actuator_1	133	50
temperature_wire_actuator_1_134	dummy-meas	temperature	wire_actuator_1	134	51
temperature_wire_actuator_1_135	dummy-meas	temperature	wire_actuator_1	135	52
temperature_wire_actuator_1_136	dummy-meas	temperature	wire_actuator_1	136	53
temperature_wire_actuator_1_137	dummy-meas	temperature	wire_actuator_1	137	60
temperature_wire_actuator_1_138	dummy-meas	temperature	wire_actuator_1	138	65
temperature_wire_actuator_1_139	dummy-meas	temperature	wire_actuator_1	139	70
hysteresis_wire_actuator_1_123	dummy-meas	hysteresis	wire_actuator_1	123	2
hysteresis_wire_actuator_1_122	dummy-meas	hysteresis	wire_actuator_1	122	2
hysteresis_wire_actuator_1_121	dummy-meas	hysteresis	wire_actuator_1	121	2
hysteresis_wire_actuator_1_120	dummy-meas	hysteresis	wire_actuator_1	120	2
hysteresis_wire_actuator_1_119	dummy-meas	hysteresis	wire_actuator_1	119	2
hysteresis_wire_actuator_1_118	dummy-meas	hysteresis	wire_actuator_1	118	2
hysteresis_wire_actuator_1_117	dummy-meas	hysteresis	wire_actuator_1	117	2
hysteresis_wire_actuator_1_116	dummy-meas	hysteresis	wire_actuator_1	116	2
hysteresis_wire_actuator_1_115	dummy-meas	hysteresis	wire_actuator_1	115	2
hysteresis_wire_actuator_1_111	dummy-meas	hysteresis	wire_actuator_1	111	2
hysteresis_wire_actuator_1_114	dummy-meas	hysteresis	wire_actuator_1	114	2
hysteresis_wire_actuator_1_113	dummy-meas	hysteresis	wire_actuator_1	113	2
hysteresis_wire_actuator_1_112	dummy-meas	hysteresis	wire_actuator_1	112	2
hysteresis_wire_actuator_1_131	dummy-meas	hysteresis	wire_actuator_1	131	1
hysteresis_wire_actuator_1_140	dummy-meas	hysteresis	wire_actuator_1	140	1
hysteresis_wire_actuator_1_141	dummy-meas	hysteresis	wire_actuator_1	141	1
hysteresis_wire_actuator_1_142	dummy-meas	hysteresis	wire_actuator_1	142	1
hysteresis_wire_actuator_1_143	dummy-meas	hysteresis	wire_actuator_1	143	1
hysteresis_wire_actuator_1_132	dummy-meas	hysteresis	wire_actuator_1	132	1
hysteresis_wire_actuator_1_133	dummy-meas	hysteresis	wire_actuator_1	133	1
hysteresis_wire_actuator_1_134	dummy-meas	hysteresis	wire_actuator_1	134	1
hysteresis_wire_actuator_1_135	dummy-meas	hysteresis	wire_actuator_1	135	1
hysteresis_wire_actuator_1_136	dummy-meas	hysteresis	wire_actuator_1	136	1
hysteresis_wire_actuator_1_137	dummy-meas	hysteresis	wire_actuator_1	137	1
hysteresis_wire_actuator_1_138	dummy-meas	hysteresis	wire_actuator_1	138	1
hysteresis_wire_actuator_1_139	dummy-meas	hysteresis	wire_actuator_1	139	1
austenite_start_temperature_wire_actuator_1_0	dummy-meas	austenite_start_temperature	wire_actuator_1	0	69.4
austenite_finish_temperature_wire_actuator_1_0	dummy-meas	austenite_finish_temperature	wire_actuator_1	0	73.9
martensite_finish_temperature_wire_actuator_1_0	dummy-meas	martensite_finish_temperature	wire_actuator_1	0	47.2
martensite_start_temperature_wire_actuator_1_0	dummy-meas	martensite_start_temperature	wire_actuator_1	0	52.4
strain_wire_actuator_NiTi#6_110	dummy-meas	strain	wire_actuator_NiTi#6	110	0
strain_wire_actuator_NiTi#6_111	dummy-meas	strain	wire_actuator_NiTi#6	111	0.005
strain_wire_actuator_NiTi#6_112	dummy-meas	strain	wire_actuator_NiTi#6	112	0.01
strain_wire_actuator_NiTi#6_113	dummy-meas	strain	wire_actuator_NiTi#6	113	0.015
strain_wire_actuator_NiTi#6_114	dummy-meas	strain	wire_actuator_NiTi#6	114	0.02
strain_wire_actuator_NiTi#6_115	dummy-meas	strain	wire_actuator_NiTi#6	115	0.025
strain_wire_actuator_NiTi#6_116	dummy-meas	strain	wire_actuator_NiTi#6	116	0.03
strain_wire_actuator_NiTi#6_117	dummy-meas	strain	wire_actuator_NiTi#6	117	0.035
strain_wire_actuator_NiTi#6_118	dummy-meas	strain	wire_actuator_NiTi#6	118	0.04
strain_wire_actuator_NiTi#6_119	dummy-meas	strain	wire_actuator_NiTi#6	119	0.045
strain_wire_actuator_NiTi#6_120	dummy-meas	strain	wire_actuator_NiTi#6	120	0.05
strain_wire_actuator_NiTi#6_121	dummy-meas	strain	wire_actuator_NiTi#6	121	0.055
strain_wire_actuator_NiTi#6_122	dummy-meas	strain	wire_actuator_NiTi#6	122	0.06
strain_wire_actuator_NiTi#6_210	dummy-meas	strain	wire_actuator_NiTi#6	210	0
strain_wire_actuator_NiTi#6_211	dummy-meas	strain	wire_actuator_NiTi#6	211	0.005
strain_wire_actuator_NiTi#6_212	dummy-meas	strain	wire_actuator_NiTi#6	212	0.01
strain_wire_actuator_NiTi#6_213	dummy-meas	strain	wire_actuator_NiTi#6	213	0.015
strain_wire_actuator_NiTi#6_214	dummy-meas	strain	wire_actuator_NiTi#6	214	0.02
strain_wire_actuator_NiTi#6_215	dummy-meas	strain	wire_actuator_NiTi#6	215	0.025
strain_wire_actuator_NiTi#6_216	dummy-meas	strain	wire_actuator_NiTi#6	216	0.03
strain_wire_actuator_NiTi#6_217	dummy-meas	strain	wire_actuator_NiTi#6	217	0.035
strain_wire_actuator_NiTi#6_218	dummy-meas	strain	wire_actuator_NiTi#6	218	0.04
strain_wire_actuator_NiTi#6_219	dummy-meas	strain	wire_actuator_NiTi#6	219	0.045
strain_wire_actuator_NiTi#6_220	dummy-meas	strain	wire_actuator_NiTi#6	220	0.05
strain_wire_actuator_NiTi#6_221	dummy-meas	strain	wire_actuator_NiTi#6	221	0.055
strain_wire_actuator_NiTi#6_222	dummy-meas	strain	wire_actuator_NiTi#6	222	0.06
strain_wire_actuator_NiTi#6_310	dummy-meas	strain	wire_actuator_NiTi#6	310	0
strain_wire_actuator_NiTi#6_311	dummy-meas	strain	wire_actuator_NiTi#6	311	0.005
strain_wire_actuator_NiTi#6_312	dummy-meas	strain	wire_actuator_NiTi#6	312	0.01
strain_wire_actuator_NiTi#6_313	dummy-meas	strain	wire_actuator_NiTi#6	313	0.015
strain_wire_actuator_NiTi#6_314	dummy-meas	strain	wire_actuator_NiTi#6	314	0.02
strain_wire_actuator_NiTi#6_315	dummy-meas	strain	wire_actuator_NiTi#6	315	0.025
strain_wire_actuator_NiTi#6_316	dummy-meas	strain	wire_actuator_NiTi#6	316	0.03
strain_wire_actuator_NiTi#6_317	dummy-meas	strain	wire_actuator_NiTi#6	317	0.035
strain_wire_actuator_NiTi#6_318	dummy-meas	strain	wire_actuator_NiTi#6	318	0.04
strain_wire_actuator_NiTi#6_319	dummy-meas	strain	wire_actuator_NiTi#6	319	0.045
strain_wire_actuator_NiTi#6_320	dummy-meas	strain	wire_actuator_NiTi#6	320	0.05
strain_wire_actuator_NiTi#6_321	dummy-meas	strain	wire_actuator_NiTi#6	321	0.055
strain_wire_actuator_NiTi#6_322	dummy-meas	strain	wire_actuator_NiTi#6	322	0.06
martensite_content_wire_actuator_NiTi#6_110	dummy-meas	martensite_content	wire_actuator_NiTi#6	110	0
martensite_content_wire_actuator_NiTi#6_111	dummy-meas	martensite_content	wire_actuator_NiTi#6	111	0
martensite_content_wire_actuator_NiTi#6_112	dummy-meas	martensite_content	wire_actuator_NiTi#6	112	0
martensite_content_wire_actuator_NiTi#6_113	dummy-meas	martensite_content	wire_actuator_NiTi#6	113	0
martensite_content_wire_actuator_NiTi#6_114	dummy-meas	martensite_content	wire_actuator_NiTi#6	114	0
martensite_content_wire_actuator_NiTi#6_115	dummy-meas	martensite_content	wire_actuator_NiTi#6	115	0
martensite_content_wire_actuator_NiTi#6_116	dummy-meas	martensite_content	wire_actuator_NiTi#6	116	0
martensite_content_wire_actuator_NiTi#6_117	dummy-meas	martensite_content	wire_actuator_NiTi#6	117	0
martensite_content_wire_actuator_NiTi#6_118	dummy-meas	martensite_content	wire_actuator_NiTi#6	118	0
martensite_content_wire_actuator_NiTi#6_119	dummy-meas	martensite_content	wire_actuator_NiTi#6	119	0
martensite_content_wire_actuator_NiTi#6_120	dummy-meas	martensite_content	wire_actuator_NiTi#6	120	0
martensite_content_wire_actuator_NiTi#6_121	dummy-meas	martensite_content	wire_actuator_NiTi#6	121	0
martensite_content_wire_actuator_NiTi#6_122	dummy-meas	martensite_content	wire_actuator_NiTi#6	122	0
martensite_content_wire_actuator_NiTi#6_210	dummy-meas	martensite_content	wire_actuator_NiTi#6	210	0.5
martensite_content_wire_actuator_NiTi#6_211	dummy-meas	martensite_content	wire_actuator_NiTi#6	211	0.5
martensite_content_wire_actuator_NiTi#6_212	dummy-meas	martensite_content	wire_actuator_NiTi#6	212	0.5
martensite_content_wire_actuator_NiTi#6_213	dummy-meas	martensite_content	wire_actuator_NiTi#6	213	0.5
martensite_content_wire_actuator_NiTi#6_214	dummy-meas	martensite_content	wire_actuator_NiTi#6	214	0.5
martensite_content_wire_actuator_NiTi#6_215	dummy-meas	martensite_content	wire_actuator_NiTi#6	215	0.5
martensite_content_wire_actuator_NiTi#6_216	dummy-meas	martensite_content	wire_actuator_NiTi#6	216	0.5
martensite_content_wire_actuator_NiTi#6_217	dummy-meas	martensite_content	wire_actuator_NiTi#6	217	0.5
martensite_content_wire_actuator_NiTi#6_218	dummy-meas	martensite_content	wire_actuator_NiTi#6	218	0.5
martensite_content_wire_actuator_NiTi#6_219	dummy-meas	martensite_content	wire_actuator_NiTi#6	219	0.5
martensite_content_wire_actuator_NiTi#6_220	dummy-meas	martensite_content	wire_actuator_NiTi#6	220	0.5
martensite_content_wire_actuator_NiTi#6_221	dummy-meas	martensite_content	wire_actuator_NiTi#6	221	0.5
martensite_content_wire_actuator_NiTi#6_222	dummy-meas	martensite_content	wire_actuator_NiTi#6	222	0.5
martensite_content_wire_actuator_NiTi#6_310	dummy-meas	martensite_content	wire_actuator_NiTi#6	310	1
martensite_content_wire_actuator_NiTi#6_311	dummy-meas	martensite_content	wire_actuator_NiTi#6	311	1
martensite_content_wire_actuator_NiTi#6_312	dummy-meas	martensite_content	wire_actuator_NiTi#6	312	1
martensite_content_wire_actuator_NiTi#6_313	dummy-meas	martensite_content	wire_actuator_NiTi#6	313	1
martensite_content_wire_actuator_NiTi#6_314	dummy-meas	martensite_content	wire_actuator_NiTi#6	314	1
martensite_content_wire_actuator_NiTi#6_315	dummy-meas	martensite_content	wire_actuator_NiTi#6	315	1
martensite_content_wire_actuator_NiTi#6_316	dummy-meas	martensite_content	wire_actuator_NiTi#6	316	1
martensite_content_wire_actuator_NiTi#6_317	dummy-meas	martensite_content	wire_actuator_NiTi#6	317	1
martensite_content_wire_actuator_NiTi#6_318	dummy-meas	martensite_content	wire_actuator_NiTi#6	318	1
martensite_content_wire_actuator_NiTi#6_319	dummy-meas	martensite_content	wire_actuator_NiTi#6	319	1
martensite_content_wire_actuator_NiTi#6_320	dummy-meas	martensite_content	wire_actuator_NiTi#6	320	1
martensite_content_wire_actuator_NiTi#6_321	dummy-meas	martensite_content	wire_actuator_NiTi#6	321	1
martensite_content_wire_actuator_NiTi#6_322	dummy-meas	martensite_content	wire_actuator_NiTi#6	322	1
martensite_content_wire_actuator_NiTi#6_500	dummy-meas	martensite_content	wire_actuator_NiTi#6	500	1
E-modulus_twinned_martensite_wire_actuator_NiTi#6_0	dummy-meas	E-modulus_twinned_martensite	wire_actuator_NiTi#6	0	22800000000
hardening_parameter_wire_actuator_NiTi#6_0	dummy-meas	hardening_parameter	wire_actuator_NiTi#6	0	150000000
plateau_finish_strain_wire_actuator_NiTi#6_0	dummy-meas	plateau_finish_strain	wire_actuator_NiTi#6	0	0.0581
plateau_start_strain_wire_actuator_NiTi#6_0	dummy-meas	plateau_start_strain	wire_actuator_NiTi#6	0	0.0088
E-modulus_de-twinned_martensite_wire_actuator_NiTi#6_0	dummy-meas	E-modulus_de-twinned_martensite	wire_actuator_NiTi#6	0	24830000000
E-modulus_austenite_wire_actuator_NiTi#6_0	dummy-meas	E-modulus_austenite	wire_actuator_NiTi#6	0	53510000000
max_blocking_stress_NiTi	dummy-meas	maximal_blocking_stress	NiTi	0	500000000
M420_1_5	dummy-meas	A00	M420	0	4.5e-11
M420_1_6	dummy-meas	B00	M420	0	-1.6e-10
M420_1_7	dummy-meas	B00	M420	0	3.55e-10
M420_1_8	dummy-meas	B00	M420	0	5.25e-10
M420_2_1	dummy-meas	B01	M420	0	-0.011
M420_2_2	dummy-meas	B01	M420	0	0.025
M420_2_3	dummy-meas	B01	M420	0	0.037
dielectric_strength_M420	dummy-meas	dielectric_strength	M420	0	2000
lattice_constant_a_Ni2MnGa_sample_1	dummy-meas	lattice_constant_a	Ni2MnGa_sample_1	0	0.59
lattice_constant_c_Ni2MnGa_sample_1	dummy-meas	lattice_constant_c	Ni2MnGa_sample_1	0	0.55
M420_1_9	dummy-meas	C00	M420	0	1.4166399999999998e-08
M420_1_10	dummy-meas	C00	M420	0	1.4166399999999998e-08
\.


--
-- Data for Name: concr_param_matrix; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.concr_param_matrix (concr_param_id, index_id, value) FROM stdin;
M420_1_1	1	1
M420_1_1	2	1
M420_1_2	1	3
M420_1_2	2	3
M420_1_3	1	1
M420_1_3	2	3
M420_1_4	1	1
M420_1_4	2	2
M420_1_5	1	5
M420_1_5	2	5
M420_1_6	1	3
M420_1_6	2	1
M420_1_7	1	3
M420_1_7	2	3
M420_1_8	1	1
M420_1_8	2	5
M420_1_9	1	1
M420_1_9	2	1
M420_1_10	1	3
M420_1_10	2	3
M420_2_1	1	3
M420_2_1	2	1
M420_2_2	1	3
M420_2_2	2	3
M420_2_3	1	1
M420_2_3	2	5
\.


--
-- Data for Name: concr_param_matrix_additional; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.concr_param_matrix_additional (concr_param_id, index_id, value, callversion) FROM stdin;
cp_mod_A11_matrix_model_1_M420_0_2	2	2	\N
cp_mod_A11_matrix_model_1_M420_0_3	2	3	\N
cp_mod_A11_matrix_model_1_M420_0_4	1	3	\N
cp_mod_A11_matrix_model_1_M420_0_4	2	3	\N
cp_mod_A11_matrix_model_1_M420_0_5	1	5	\N
cp_mod_A11_matrix_model_1_M420_0_5	2	5	\N
cp_mod_B00_matrix_model_2_M420_0_1	1	3	\N
cp_mod_B00_matrix_model_2_M420_0_1	2	1	\N
cp_mod_B00_matrix_model_2_M420_0_2	1	3	\N
cp_mod_B00_matrix_model_2_M420_0_2	2	3	\N
cp_mod_B00_matrix_model_2_M420_0_3	1	1	\N
cp_mod_B00_matrix_model_2_M420_0_3	2	5	\N
cp_mod_B11_matrix_model_2_M420_0_1	1	3	\N
cp_mod_B11_matrix_model_2_M420_0_1	2	1	\N
cp_mod_B11_matrix_model_2_M420_0_2	1	3	\N
cp_mod_B11_matrix_model_2_M420_0_2	2	3	\N
cp_mod_B11_matrix_model_2_M420_0_3	1	1	\N
cp_mod_B11_matrix_model_2_M420_0_3	2	5	\N
cp_mod_C00_matrix_model_3_M420_0_1	1	1	\N
cp_mod_C00_matrix_model_3_M420_0_1	2	1	\N
cp_mod_C00_matrix_model_3_M420_0_2	1	3	\N
cp_mod_C00_matrix_model_3_M420_0_2	2	3	\N
cp_mod_C11_matrix_model_3_M420_0_1	1	1	\N
cp_mod_C11_matrix_model_3_M420_0_1	2	1	\N
cp_mod_C11_matrix_model_3_M420_0_2	1	3	\N
cp_mod_C11_matrix_model_3_M420_0_2	2	3	\N
cp_mod_A00_matrix_model_4_M420_0_1	1	1	\N
cp_mod_A00_matrix_model_4_M420_0_1	2	1	\N
cp_mod_A00_matrix_model_4_M420_0_2	1	1	\N
cp_mod_A00_matrix_model_4_M420_0_2	2	2	\N
cp_mod_A00_matrix_model_4_M420_0_3	1	1	\N
cp_mod_A00_matrix_model_4_M420_0_3	2	3	\N
cp_mod_A00_matrix_model_4_M420_0_4	1	3	\N
cp_mod_A00_matrix_model_4_M420_0_4	2	3	\N
cp_mod_A00_matrix_model_4_M420_0_5	1	5	\N
cp_mod_A00_matrix_model_4_M420_0_5	2	5	\N
cp_mod_A11_matrix_model_4_M420_0_1	1	1	\N
cp_mod_A11_matrix_model_4_M420_0_1	2	1	\N
cp_mod_A11_matrix_model_4_M420_0_2	1	1	\N
cp_mod_A11_matrix_model_4_M420_0_2	2	2	\N
cp_mod_A11_matrix_model_4_M420_0_3	1	1	\N
cp_mod_A11_matrix_model_4_M420_0_3	2	3	\N
cp_mod_A11_matrix_model_4_M420_0_4	1	3	\N
cp_mod_A11_matrix_model_4_M420_0_4	2	3	\N
cp_mod_A11_matrix_model_4_M420_0_5	1	5	\N
cp_mod_A11_matrix_model_4_M420_0_5	2	5	\N
cp_mod_B00_matrix_model_5_M420_0_1	1	3	\N
cp_mod_B00_matrix_model_5_M420_0_1	2	1	\N
cp_mod_B00_matrix_model_5_M420_0_2	1	3	\N
cp_mod_B00_matrix_model_5_M420_0_2	2	3	\N
cp_mod_B00_matrix_model_5_M420_0_3	1	1	\N
cp_mod_B00_matrix_model_5_M420_0_3	2	5	\N
cp_mod_B11_matrix_model_5_M420_0_1	1	3	\N
cp_mod_B11_matrix_model_5_M420_0_1	2	1	\N
cp_mod_B11_matrix_model_5_M420_0_2	1	3	\N
cp_mod_B11_matrix_model_5_M420_0_2	2	3	\N
cp_mod_B11_matrix_model_5_M420_0_3	1	1	\N
cp_mod_B11_matrix_model_5_M420_0_3	2	5	\N
cp_mod_C00_matrix_model_6_M420_0_1	1	1	\N
cp_mod_C00_matrix_model_6_M420_0_1	2	1	\N
cp_mod_C00_matrix_model_6_M420_0_2	1	3	\N
cp_mod_C00_matrix_model_6_M420_0_2	2	3	\N
cp_mod_C11_matrix_model_6_M420_0_1	1	1	\N
cp_mod_C11_matrix_model_6_M420_0_1	2	1	\N
cp_mod_C11_matrix_model_6_M420_0_2	1	3	\N
cp_mod_C11_matrix_model_6_M420_0_2	2	3	\N
cp_mod_A01_matrix_model_1_M420_0_1	1	1	\N
cp_mod_A01_matrix_model_1_M420_0_1	2	1	\N
cp_mod_A01_matrix_model_1_M420_0_2	1	1	\N
cp_mod_A01_matrix_model_1_M420_0_2	2	2	\N
cp_mod_A01_matrix_model_1_M420_0_3	1	1	\N
cp_mod_A01_matrix_model_1_M420_0_3	2	3	\N
cp_mod_A01_matrix_model_1_M420_0_4	1	3	\N
cp_mod_A01_matrix_model_1_M420_0_4	2	3	\N
cp_mod_A01_matrix_model_1_M420_0_5	1	5	\N
cp_mod_A01_matrix_model_1_M420_0_5	2	5	\N
cp_mod_B01_matrix_model_2_M420_0_1	1	3	\N
cp_mod_B01_matrix_model_2_M420_0_1	2	1	\N
cp_mod_B01_matrix_model_2_M420_0_2	1	3	\N
cp_mod_B01_matrix_model_2_M420_0_2	2	3	\N
cp_mod_B01_matrix_model_2_M420_0_3	1	1	\N
cp_mod_B01_matrix_model_2_M420_0_3	2	5	\N
cp_mod_C01_matrix_model_3_M420_0_1	1	1	\N
cp_mod_C01_matrix_model_3_M420_0_1	2	1	\N
cp_mod_C01_matrix_model_3_M420_0_2	1	3	\N
cp_mod_C01_matrix_model_3_M420_0_2	2	3	\N
cp_mod_A10_matrix_model_4_M420_0_1	1	1	\N
cp_mod_A10_matrix_model_4_M420_0_1	2	1	\N
cp_mod_A10_matrix_model_4_M420_0_2	1	1	\N
cp_mod_A10_matrix_model_4_M420_0_2	2	2	\N
cp_mod_A10_matrix_model_4_M420_0_3	1	1	\N
cp_mod_A10_matrix_model_4_M420_0_3	2	3	\N
cp_mod_A10_matrix_model_4_M420_0_4	1	3	\N
cp_mod_A10_matrix_model_4_M420_0_4	2	3	\N
cp_mod_A10_matrix_model_4_M420_0_5	1	5	\N
cp_mod_A10_matrix_model_4_M420_0_5	2	5	\N
cp_mod_B10_matrix_model_5_M420_0_1	1	3	\N
cp_mod_B10_matrix_model_5_M420_0_1	2	1	\N
cp_mod_B10_matrix_model_5_M420_0_2	1	3	\N
cp_mod_B10_matrix_model_5_M420_0_2	2	3	\N
cp_mod_B10_matrix_model_5_M420_0_3	1	1	\N
cp_mod_B10_matrix_model_5_M420_0_3	2	5	\N
cp_mod_C10_matrix_model_6_M420_0_1	1	1	\N
cp_mod_C10_matrix_model_6_M420_0_1	2	1	\N
cp_mod_C10_matrix_model_6_M420_0_2	1	3	\N
cp_mod_C10_matrix_model_6_M420_0_2	2	3	\N
cp_mod_A10_matrix_model_1_M420_0_1	1	1	\N
cp_mod_A10_matrix_model_1_M420_0_1	2	1	\N
cp_mod_A10_matrix_model_1_M420_0_2	1	1	\N
cp_mod_A10_matrix_model_1_M420_0_2	2	2	\N
cp_mod_A10_matrix_model_1_M420_0_3	1	1	\N
cp_mod_A10_matrix_model_1_M420_0_3	2	3	\N
cp_mod_A10_matrix_model_1_M420_0_4	1	3	\N
cp_mod_A10_matrix_model_1_M420_0_4	2	3	\N
cp_mod_A10_matrix_model_1_M420_0_5	1	5	\N
cp_mod_A10_matrix_model_1_M420_0_5	2	5	\N
cp_mod_B10_matrix_model_2_M420_0_1	1	3	\N
cp_mod_B10_matrix_model_2_M420_0_1	2	1	\N
cp_mod_B10_matrix_model_2_M420_0_2	1	3	\N
cp_mod_B10_matrix_model_2_M420_0_2	2	3	\N
cp_mod_B10_matrix_model_2_M420_0_3	1	1	\N
cp_mod_B10_matrix_model_2_M420_0_3	2	5	\N
cp_mod_C10_matrix_model_3_M420_0_1	1	1	\N
cp_mod_C10_matrix_model_3_M420_0_1	2	1	\N
cp_mod_C10_matrix_model_3_M420_0_2	1	3	\N
cp_mod_C10_matrix_model_3_M420_0_2	2	3	\N
cp_mod_A01_matrix_model_4_M420_0_1	1	1	\N
cp_mod_A01_matrix_model_4_M420_0_1	2	1	\N
cp_mod_A01_matrix_model_4_M420_0_2	1	1	\N
cp_mod_A01_matrix_model_4_M420_0_2	2	2	\N
cp_mod_A01_matrix_model_4_M420_0_3	1	1	\N
cp_mod_A01_matrix_model_4_M420_0_3	2	3	\N
cp_mod_A01_matrix_model_4_M420_0_4	1	3	\N
cp_mod_A01_matrix_model_4_M420_0_4	2	3	\N
cp_mod_A01_matrix_model_4_M420_0_5	1	5	\N
cp_mod_A01_matrix_model_4_M420_0_5	2	5	\N
cp_mod_B01_matrix_model_5_M420_0_1	1	3	\N
cp_mod_B01_matrix_model_5_M420_0_1	2	1	\N
cp_mod_B01_matrix_model_5_M420_0_2	1	3	\N
cp_mod_B01_matrix_model_5_M420_0_2	2	3	\N
cp_mod_B01_matrix_model_5_M420_0_3	1	1	\N
cp_mod_B01_matrix_model_5_M420_0_3	2	5	\N
cp_mod_C01_matrix_model_6_M420_0_1	1	1	\N
cp_mod_C01_matrix_model_6_M420_0_1	2	1	\N
cp_mod_C01_matrix_model_6_M420_0_2	1	3	\N
cp_mod_C01_matrix_model_6_M420_0_2	2	3	\N
cp_mod_A00_matrix_model_1_M420_0_1	1	1	\N
cp_mod_A00_matrix_model_1_M420_0_1	2	1	\N
cp_mod_A00_matrix_model_1_M420_0_2	1	1	\N
cp_mod_A00_matrix_model_1_M420_0_2	2	2	\N
cp_mod_A00_matrix_model_1_M420_0_3	1	1	\N
cp_mod_A00_matrix_model_1_M420_0_3	2	3	\N
cp_mod_A00_matrix_model_1_M420_0_4	1	3	\N
cp_mod_A00_matrix_model_1_M420_0_4	2	3	\N
cp_mod_A00_matrix_model_1_M420_0_5	1	5	\N
cp_mod_A00_matrix_model_1_M420_0_5	2	5	\N
cp_mod_A11_matrix_model_1_M420_0_1	1	1	\N
cp_mod_A11_matrix_model_1_M420_0_1	2	1	\N
cp_mod_A11_matrix_model_1_M420_0_2	1	1	\N
cp_mod_A11_matrix_model_1_M420_0_3	1	1	\N
\.


--
-- Data for Name: constants; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.constants (param_id, value) FROM stdin;
vacuum_permittivity	0.000000000008854
vacuum_permeability	0.0000012566
\.


--
-- Data for Name: curve_condition; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.curve_condition (curve_id, param_id, value, in_number, operator) FROM stdin;
magnetic_curve_hard_axis	easy_hard_axis	1	-1	\N
magnetic_curve_easy_axis	easy_hard_axis	0	-1	\N
hysteresis_martensite_content	hysteresis	\N	-100	\N
stress-strain_curve_martensite	martensite_content	\N	-100	\N
\.


--
-- Data for Name: curve_params; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.curve_params (curve_id, param_id) FROM stdin;
magnetic_curve_easy_axis	magnetic_field_strength
magnetic_curve_easy_axis	magnetization
magnetic_curve_easy_axis	flux_density
magnetic_curve_easy_axis	magnetic_polarization
magnetic_curve_easy_axis	relative_permeability
magnetic_curve_hard_axis	magnetic_field_strength
magnetic_curve_hard_axis	magnetization
magnetic_curve_hard_axis	flux_density
magnetic_curve_hard_axis	magnetic_polarization
magnetic_curve_hard_axis	relative_permeability
hysteresis_martensite_content	temperature
hysteresis_martensite_content	martensite_content
stress-strain_curve_martensite	mechanic_stress
stress-strain_curve_martensite	strain
stress-strain_curve_elastomer	mechanic_stress
stress-strain_curve_elastomer	strain
curve_magnetic_stress	magnetic_stress
curve_magnetic_stress	magnetic_field_strength
\.


--
-- Data for Name: material; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.material (mat_id, mat_name, mat_type) FROM stdin;
Ni2MnGa	Ni2MnGa	MSMA
NiTi	NiTi	SMA
M420	M420	PC
elastomer_1	Wacker_Elastosil2030	DE
\.


--
-- Data for Name: material_condition; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.material_condition (mod_id, mat_type, abstract_object_id) FROM stdin;
calc_young_modulus_initial_slope	DE	\N
calc_neo_hookean_initial_slope	DE	\N
calc_young_modulus_minR	DE	\N
calc_young_modulus_Neo_Hookean_minR	DE	\N
calc_yeoh_minR	DE	\N
\.


--
-- Data for Name: material_types; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.material_types (mat_type) FROM stdin;
DE
SMA
MSMA
PC
\.


--
-- Data for Name: measinput; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.measinput (meas_id, param_id) FROM stdin;
\.


--
-- Data for Name: measoutput; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.measoutput (meas_id, param_id) FROM stdin;
\.


--
-- Data for Name: measurement; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.measurement (meas_id, meas_name) FROM stdin;
meas_flux_over_field	measurement of flux density over magnetic field
tensile_test	tensile test
base_value	 
measured_at	 given data
dummy-data	 given data
\.


--
-- Data for Name: model; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.model (mod_id, mod_name, mod_type, mod_equation, function_name) FROM stdin;
electrostatic_pressure_model	electrostatic pressure model	simple_calculation	σ_el = ϵ_r*ϵ_0*E^2	mod_el_pres
calc_magnetization	calc magnetization based on characteristic curve	simple_calculation	M = B μ_0 -H	calc_magnetization
calc_flux_density	calc flux density based on characteristic curve	simple_calculation	B = μ_0 μ_r H	calc_flux_density
calc_mag_polarization	calc magnetic polarization based on characteristic curve	simple_calculation	J = M μ_0	calc_mag_polarization
calc_rel_perm	calc relative permeability based on characteristic curve	simple_calculation	μ_r = 1 + J/(H μ_0)	calc_rel_perm
calc_mag_energy_density	calc of the magnetic energy density	agg_integral	g_{n+1}=g_n+(H_{n+1}-H_n)*(J_{n+1}-J_n)/2	calc_mag_energy_density
calc_young_modulus_initial_slope	parameter identification with initial slope	simple_calculation	Y = σ_ini / ε_ini	ym_init_slope
calc_neo_hookean_initial_slope	parameter identifikation Neo Hookean with initial slope	simple_calculation	Y_NH = 3σ_i / (ε_i^2 - 1/ε_i)	calc_neo_hookean_initial_slope
calc_young_modulus_minR	calc of young modulus with error minimization	aggregate	minR over Y = σ / ε	calc_young_modulus_minr
calc_young_modulus_Neo_Hookean_minR	calc of young modulus Neo Hookean with error minimization	aggregate	minR over Y_NH = 3σ_i / (λ_i^2 - 1/λ_i)	calc_neo_hook_minr
calc_yeoh_minR	calc of the Yeoh parameter vector with error minimization	aggregate	\N	calc_yeoh_minr
calc_martensite_asc	martensite content in relation to temperature, ascending	simple_calculation	1/2* cos(PI*(T-Mf)/(Ms-Mf))+1/2	calc_martensite_asc
calc_martensite_desc	martensite content in relation to temperature, descending	simple_calculation	1/2* cos(PI*(T-As)/(Af-As))+1/2	calc_martensite_desc
1d_non_linear_mech_mod_v1	one-dimensional non-linear mechanical model, 0<ε<ε_ps	simple_calculation	σ = ε⋅(E_A - (E_A -E_M) * ξ)	mech_mod_v1
1d_non_linear_mech_mod_v2	one-dimensional non-linear mechanical model, ε_ps<ε<ε_pf	simple_calculation	σ = ε_ps⋅(E_A - (E_A -E_M) * ξ)+(ε-ε_ps)*H	mech_mod_v2
1d_non_linear_mech_mod_v3	one-dimensional non-linear mechanical model, ε_pf<ε	simple_calculation	σ = ε_ps⋅(E_A - (E_A -E_M) * ξ)+(ε_pf-ε_ps)*H+(ε-ε_pf)*E_D	mech_mod_v3
calc_magnetic_stress	calc of magnetic stress with magnetic curve	complex_sum	σ_M = (ge-gh)/ϵ_0	calc_magn_stress
calc_max_magn_stress	max magnetic stress as maximum of characteristic curve	aggregate_simple	σ_M = max(σ)	max
linear_interpolation_initial_stress	linear interpolation	linear_interpolation	y=(y_1-y_2)/(x_1-x_2)*(x-x_1) + y_1)	linear_interpol_initial_stress
calc_max_block_load_hold	calc of maximum blocking stress hold and load	simple_calculation	(1.) σ_B,max,load = σ_mag,max - σ_tw (2.) σ_B,max,hold = σ_mag,max + σ_tw	calc_max_block_stress_lh
calc_max_blocking_stress_PC	calc maximum blocking stress for PC materials	simple_calculation	σ_B-max = d_33,1/s^E_33,1*E_BFS	calc_max_block_pc
calc_maximum_strain_mag	calc maximum strain in the magnetic field via lattice constants	simple_calculation	ϵ_0 = 1-c/a	calc_max_strain_mag
matrix_model_2	Matrixtransformation	matrixoperation	 B^{\\bar m e} = B^{me}(A^{me})^{-1}	matrix_func_2
matrix_model_1	Matrixtransformation	matrixoperation	A^{\\bar m e} = (A^{me})^{-1}	matrix_func_1
matrix_model_3	Matrixtransformation	matrixoperation	C^{\\bar m e} = C^{me} - [e][m] B^{me} (A^{me})^{-1}(B^{me})^T	matrix_func_3
matrix_model_4	Matrixtransformation	matrixoperation	A^{m \\bar e} = A^{me} - [e][m] (B^{me})^T (C^{me})^{-1} B^{me}	matrix_func_4
matrix_model_5	Matrixtransformation	matrixoperation	B^{m \\bar e} = (C^{me})^{-1} B^{me}	matrix_func_5
matrix_model_6	Matrixtransformation	matrixoperation	C^{m \\bar e} &= (C^{me})^{-1}	matrix_func_6
\.


--
-- Data for Name: model_condition; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.model_condition (mod_id, param_id, value, in_number, operator) FROM stdin;
calc_martensite_desc	hysteresis	2	-1	=
calc_martensite_asc	hysteresis	1	-1	=
\.


--
-- Data for Name: model_condition_compare; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.model_condition_compare (mod_id, param_id, in_number, operator, compare_position) FROM stdin;
1d_non_linear_mech_mod_v1	strain	-1	<	1
1d_non_linear_mech_mod_v1	plateau_start_strain	-1	<	2
1d_non_linear_mech_mod_v2	strain	-1	>	1
1d_non_linear_mech_mod_v2	plateau_start_strain	-1	>	2
1d_non_linear_mech_mod_v2	strain	-2	<	1
1d_non_linear_mech_mod_v2	plateau_finish_strain	-2	<	2
1d_non_linear_mech_mod_v3	plateau_finish_strain	-1	<	1
1d_non_linear_mech_mod_v3	strain	-1	<	2
\.


--
-- Data for Name: model_condition_interpol; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.model_condition_interpol (mod_id, param_id, value, in_number) FROM stdin;
\.


--
-- Data for Name: modelinput; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.modelinput (mod_id, param_id, in_number) FROM stdin;
electrostatic_pressure_model	relative_permittivity	1
electrostatic_pressure_model	vacuum_permittivity	2
electrostatic_pressure_model	dielectric_strength	3
calc_max_block_load_hold	max_magnetic_stress	1
calc_max_block_load_hold	twinning_stress	2
calc_magnetization	magnetic_field_strength	1
calc_magnetization	flux_density	2
calc_magnetization	vacuum_permeability	3
calc_flux_density	vacuum_permeability	1
calc_flux_density	relative_permeability	2
calc_flux_density	magnetic_field_strength	3
calc_mag_polarization	magnetization	1
calc_mag_polarization	vacuum_permeability	2
calc_rel_perm	magnetic_polarization	1
calc_rel_perm	magnetic_field_strength	2
calc_rel_perm	vacuum_permeability	3
calc_magnetic_stress	magnetic_field_strength	1
calc_magnetic_stress	easy_hard_axis	3
calc_magnetic_stress	maximum_strain_mag	4
calc_magnetic_stress	mag_energy_density	2
calc_mag_energy_density	magnetic_field_strength	1
calc_mag_energy_density	easy_hard_axis	3
calc_mag_energy_density	magnetic_polarization	2
calc_young_modulus_minR	mechanic_stress	1
calc_young_modulus_minR	strain	2
calc_young_modulus_Neo_Hookean_minR	strain	2
calc_young_modulus_Neo_Hookean_minR	mechanic_stress	1
calc_yeoh_minR	mechanic_stress	1
calc_yeoh_minR	strain	2
matrix_model_1	A	1
matrix_model_2	A	1
matrix_model_2	B	2
matrix_model_3	A	1
matrix_model_3	B	2
matrix_model_3	C	3
matrix_model_4	A	1
matrix_model_4	B	2
matrix_model_4	C	3
matrix_model_5	B	1
matrix_model_5	C	2
matrix_model_6	C	1
calc_martensite_asc	temperature	1
calc_martensite_asc	martensite_finish_temperature	2
calc_martensite_asc	martensite_start_temperature	3
calc_martensite_desc	temperature	1
calc_martensite_desc	austenite_start_temperature	2
calc_martensite_desc	austenite_finish_temperature	3
1d_non_linear_mech_mod_v1	strain	1
1d_non_linear_mech_mod_v2	strain	1
1d_non_linear_mech_mod_v3	strain	1
1d_non_linear_mech_mod_v1	E-modulus_austenite	2
1d_non_linear_mech_mod_v2	E-modulus_austenite	2
1d_non_linear_mech_mod_v3	E-modulus_austenite	2
1d_non_linear_mech_mod_v1	E-modulus_twinned_martensite	3
1d_non_linear_mech_mod_v2	E-modulus_twinned_martensite	3
1d_non_linear_mech_mod_v3	E-modulus_twinned_martensite	3
1d_non_linear_mech_mod_v1	martensite_content	4
1d_non_linear_mech_mod_v2	martensite_content	4
1d_non_linear_mech_mod_v3	martensite_content	4
1d_non_linear_mech_mod_v2	plateau_start_strain	5
1d_non_linear_mech_mod_v3	plateau_start_strain	5
1d_non_linear_mech_mod_v2	hardening_parameter	6
1d_non_linear_mech_mod_v3	hardening_parameter	6
1d_non_linear_mech_mod_v3	plateau_finish_strain	7
1d_non_linear_mech_mod_v3	E-modulus_de-twinned_martensite	8
calc_max_magn_stress	magnetic_stress	1
calc_max_blocking_stress_PC	B00	1
calc_max_blocking_stress_PC	A00	2
calc_max_blocking_stress_PC	dielectric_strength	3
calc_neo_hookean_initial_slope	initial_strain	2
calc_neo_hookean_initial_slope	initial_stress	1
calc_young_modulus_initial_slope	initial_strain	2
calc_young_modulus_initial_slope	initial_stress	1
linear_interpolation_initial_stress	initial_strain	1
linear_interpolation_initial_stress	strain	2
linear_interpolation_initial_stress	mechanic_stress	3
calc_max_magn_stress	magnetic_field_strength	2
calc_maximum_strain_mag	lattice_constant_c	1
calc_maximum_strain_mag	lattice_constant_a	2
\.


--
-- Data for Name: modelinput_matrix; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.modelinput_matrix (mod_id, param_id, index_id, value) FROM stdin;
calc_max_blocking_stress_PC	A00	1	3
calc_max_blocking_stress_PC	A00	2	3
calc_max_blocking_stress_PC	B00	2	3
calc_max_blocking_stress_PC	B00	1	3
\.


--
-- Data for Name: modeloutput; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.modeloutput (mod_id, param_id, out_number) FROM stdin;
electrostatic_pressure_model	maximal_blocking_stress	1
calc_max_block_load_hold	maximal_blocking_stress_load	1
calc_max_block_load_hold	maximal_blocking_stress_hold	2
calc_magnetization	magnetization	1
calc_mag_polarization	magnetic_polarization	1
calc_flux_density	flux_density	1
calc_rel_perm	relative_permeability	1
calc_magnetic_stress	magnetic_stress	1
calc_mag_energy_density	mag_energy_density	1
calc_young_modulus_initial_slope	young_modulus	1
calc_neo_hookean_initial_slope	young_modulus_neo_Hookean	1
calc_young_modulus_minR	young_modulus	1
calc_young_modulus_Neo_Hookean_minR	young_modulus_neo_Hookean	1
calc_yeoh_minR	Yeoh_vector_1	1
calc_yeoh_minR	Yeoh_vector_2	2
calc_yeoh_minR	Yeoh_vector_3	3
matrix_model_1	A	1
matrix_model_2	B	1
matrix_model_3	C	1
matrix_model_4	A	1
matrix_model_5	B	1
matrix_model_6	C	1
calc_martensite_desc	martensite_content	1
calc_martensite_asc	martensite_content	1
1d_non_linear_mech_mod_v1	mechanic_stress	1
1d_non_linear_mech_mod_v2	mechanic_stress	1
1d_non_linear_mech_mod_v3	mechanic_stress	1
calc_max_magn_stress	max_magnetic_stress	1
calc_max_blocking_stress_PC	maximal_blocking_stress	1
linear_interpolation_initial_stress	initial_stress	1
calc_maximum_strain_mag	maximum_strain_mag	1
\.


--
-- Data for Name: modeloutput_pek_bits; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.modeloutput_pek_bits (mod_id, "bit") FROM stdin;
matrix_model_1	m
matrix_model_2	m
matrix_model_3	m
matrix_model_4	e
matrix_model_5	e
matrix_model_6	e
\.


--
-- Data for Name: modeltype; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.modeltype (modeltype) FROM stdin;
simple_calculation
complex_sum
agg_integral
aggregate
matrixoperation
aggregate_simple
linear_interpolation
\.


--
-- Data for Name: param; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.param (param_id, param_name, param_unit, param_symbol) FROM stdin;
electric_field	eletric field	V/m	E
vacuum_permittivity	vacuum permittivity	F/m	ε_0
relative_permittivity	relative permittivity	dimensionless quantity	ε_r
twinning_stress	twinning stress	Pa	σ_tw
magnetic_field_strength	magnetic field strength	A/m	H
magnetization	magnetization	A/m	M
flux_density	flux density	T	B
magnetic_polarization	magnetic polarization	T	J
relative_permeability	relative permeability	dimensionless quantity	μ_r
vacuum_permeability	vacuum permeability	H/m	μ_0
easy_hard_axis	easy/hard axis	dimensionless quantity	undefined
magnetic_stress	magnetic stress	Pa	σ_M
maximum_strain_mag	maximum strain in the magnetic field	dimensionless quantity	ϵ_0
mag_energy_density	free magnetization energy density	J/m^3	g
young_modulus_neo_Hookean	Young modulus according to Neo Hookean model	Pa	Y_NH
young_modulus	Young modulus	Pa	Y
Yeoh_vector_3	Yeoh vector, index 3	Pa	C_Yeoh[3]
Yeoh_vector_2	Yeoh vector, index 2	Pa	C_Yeoh[2]
Yeoh_vector_1	Yeoh vector, index 1	Pa	C_Yeoh[1]
strain	strain	dimensionless quantity	ε
mechanic_stress	mechanic stress	Pa	σ
A	\N	undefined	undefined
B	\N	undefined	undefined
C	\N	undefined	undefined
A00	mechanical compliance at constant electric field	m^2/N	s^E
C10	permittivity at constant mechanic strain	F/m	ϵ^S
C11	inverse permittivity at constant mechanic strain	m/F	β^S
C01	inverse permittivity constant mechanic stress	m/F	β^T
A01	mechanical compliance at constant flux density	m^2/N	s^D
B10	piezoelectric strain constant	N/(Vm)	e
B11	piezoelectric coupling	N/C	h
B00	piezoelectric charge constant	C/N	d
B01	piezoelectric stress constant	Vm/N	g
A10	mechanical stiffness at constant electric field	N/m^2	c^E
A11	mechanical stiffness at constant flux density	N/m^2	c^D
C00	permittivity at constant mechanic stress	F/m	ϵ^T
martensite_content	martensite content	dimensionless quantity	ξ
temperature	temperature	°C	T
martensite_start_temperature	martensite start temperature	°C	M_s
martensite_finish_temperature	martensite finish temperature	°C	M_f
austenite_start_temperature	austenite start temperature	°C	A_s
austenite_finish_temperature	austenite finish temperature	°C	A_f
hysteresis	hysteresis	dimensionless quantity	undefined
E-modulus_austenite	E-modulus austenite	Pa	E_A
E-modulus_twinned_martensite	E-modulus twinned martensite	Pa	E_M
E-modulus_de-twinned_martensite	E-modulus de-twinned martensite	Pa	E_D
hardening_parameter	hardening parameter	Pa	H
plateau_start_strain	plateau start strain	dimensionless quantity	ε_ps
plateau_finish_strain	plateau finish strain	dimensionless quantity	ε_pf
dielectric_strength	dielectric strength	V/m	E_BFS
initial_strain	initial strain	dimensionless quantity	σ_ini
initial_stress	initial stress	Pa	ε_ini
maximal_blocking_stress	maximum blocking stress	Pa	σ_B
maximal_blocking_stress_hold	maximum blocking stress hold	Pa	σ_B,max,hold
maximal_blocking_stress_load	maximum blocking stress load	Pa	σ_B,max,load
max_magnetic_stress	maximum magnetic stress	Pa	σ_M
lattice_constant_a	lattice constant a	nm	a
lattice_constant_c	lattice constant c	nm	c
\.


--
-- Data for Name: param_piezo; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.param_piezo (param_id, index_id, param_matrix, value) FROM stdin;
A00	m	A	0
A00	e	A	0
A01	m	A	0
A01	e	A	1
A10	m	A	1
A10	e	A	0
A11	m	A	1
A11	e	A	1
B00	m	B	0
B00	e	B	0
B01	m	B	0
B01	e	B	1
B10	m	B	1
B10	e	B	0
B11	m	B	1
B11	e	B	1
C00	m	C	0
C00	e	C	0
C01	m	C	0
C01	e	C	1
C10	m	C	1
C10	e	C	0
C11	m	C	1
C11	e	C	1
\.


--
-- Data for Name: param_semantic; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.param_semantic (param_id, value, semantic) FROM stdin;
easy_hard_axis	0	easy axis
easy_hard_axis	1	hard axis
hysteresis	1	ascending
hysteresis	2	descending
\.


--
-- Data for Name: sample; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.sample (sample_id, mat_id, description, sample_name) FROM stdin;
Ni2MnGa_sample_1	Ni2MnGa	\N	Ni2MnGa_sample_1
sample1_elastomer1	elastomer_1	\N	sample_Wacker_Elastosil2030
\.


--
-- Data for Name: temp_step_view; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.temp_step_view (concr_param_id, mod_meas_id, param_id, mat_sample_id, meas_time, value, mod_id) FROM stdin;
\.


--
-- Data for Name: temp_view; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.temp_view (concr_param_id, mod_meas_id, param_id, mat_sample_id, meas_time, value, mod_id, step) FROM stdin;
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_11	\N	easy_hard_axis	Ni2MnGa_sample_1	11	1	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_115	\N	easy_hard_axis	Ni2MnGa_sample_1	115	0	\N	0
strain_wire_actuator_NiTi#6_221	\N	strain	wire_actuator_NiTi#6	221	0.055	\N	0
mechanic_stress_sample1_elastomer1_tens_1_333	\N	mechanic_stress	sample1_elastomer1	333	1068300	\N	0
martensite_content_wire_actuator_NiTi#6_116	\N	martensite_content	wire_actuator_NiTi#6	116	0	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_30	\N	easy_hard_axis	Ni2MnGa_sample_1	30	1	\N	0
strain_sample1_elastomer1_tens_1_377	\N	strain	sample1_elastomer1	377	1.88	\N	0
mechanic_stress_sample1_elastomer1_tens_1_372	\N	mechanic_stress	sample1_elastomer1	372	1198300	\N	0
strain_wire_actuator_NiTi#6_314	\N	strain	wire_actuator_NiTi#6	314	0.02	\N	0
strain_sample1_elastomer1_tens_1_165	\N	strain	sample1_elastomer1	165	0.82	\N	0
strain_sample1_elastomer1_tens_1_153	\N	strain	sample1_elastomer1	153	0.76	\N	0
mechanic_stress_sample1_elastomer1_tens_1_286	\N	mechanic_stress	sample1_elastomer1	286	910300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_28	\N	mechanic_stress	sample1_elastomer1	28	147300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_130	\N	mechanic_stress	sample1_elastomer1	130	479300	\N	0
strain_sample1_elastomer1_tens_1_62	\N	strain	sample1_elastomer1	62	0.305	\N	0
strain_sample1_elastomer1_tens_1_15	\N	strain	sample1_elastomer1	15	0.07	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_134	\N	magnetic_field_strength	Ni2MnGa_sample_1	134	330000	\N	0
strain_sample1_elastomer1_tens_1_225	\N	strain	sample1_elastomer1	225	1.12	\N	0
mechanic_stress_sample1_elastomer1_tens_1_23	\N	mechanic_stress	sample1_elastomer1	23	123300	\N	0
strain_wire_actuator_NiTi#6_216	\N	strain	wire_actuator_NiTi#6	216	0.03	\N	0
mechanic_stress_sample1_elastomer1_tens_1_46	\N	mechanic_stress	sample1_elastomer1	46	224300	\N	0
strain_sample1_elastomer1_tens_1_9	\N	strain	sample1_elastomer1	9	0.04	\N	0
mechanic_stress_sample1_elastomer1_tens_1_58	\N	mechanic_stress	sample1_elastomer1	58	269300	\N	0
temperature_wire_actuator_1_139	\N	temperature	wire_actuator_1	139	70	\N	0
strain_sample1_elastomer1_tens_1_190	\N	strain	sample1_elastomer1	190	0.945	\N	0
strain_sample1_elastomer1_tens_1_327	\N	strain	sample1_elastomer1	327	1.63	\N	0
temperature_wire_actuator_1_140	\N	temperature	wire_actuator_1	140	71	\N	0
mechanic_stress_sample1_elastomer1_tens_1_52	\N	mechanic_stress	sample1_elastomer1	52	247300	\N	0
strain_sample1_elastomer1_tens_1_145	\N	strain	sample1_elastomer1	145	0.72	\N	0
strain_sample1_elastomer1_tens_1_107	\N	strain	sample1_elastomer1	107	0.53	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_162	\N	magnetic_field_strength	Ni2MnGa_sample_1	162	610000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_54	\N	easy_hard_axis	Ni2MnGa_sample_1	54	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_36	\N	magnetic_field_strength	Ni2MnGa_sample_1	36	350000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_137	\N	magnetic_field_strength	Ni2MnGa_sample_1	137	360000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_59	\N	easy_hard_axis	Ni2MnGa_sample_1	59	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_156	\N	mechanic_stress	sample1_elastomer1	156	547300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_33	\N	mechanic_stress	sample1_elastomer1	33	170300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_382	\N	mechanic_stress	sample1_elastomer1	382	1233300	\N	0
strain_sample1_elastomer1_tens_1_123	\N	strain	sample1_elastomer1	123	0.61	\N	0
mechanic_stress_sample1_elastomer1_tens_1_21	\N	mechanic_stress	sample1_elastomer1	21	113300	\N	0
strain_sample1_elastomer1_tens_1_223	\N	strain	sample1_elastomer1	223	1.11	\N	0
strain_sample1_elastomer1_tens_1_139	\N	strain	sample1_elastomer1	139	0.69	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_61	\N	easy_hard_axis	Ni2MnGa_sample_1	61	1	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_12	\N	easy_hard_axis	Ni2MnGa_sample_1	12	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_142	\N	magnetic_field_strength	Ni2MnGa_sample_1	142	410000	\N	0
strain_sample1_elastomer1_tens_1_319	\N	strain	sample1_elastomer1	319	1.59	\N	0
mechanic_stress_sample1_elastomer1_tens_1_54	\N	mechanic_stress	sample1_elastomer1	54	254300	\N	0
strain_sample1_elastomer1_tens_1_349	\N	strain	sample1_elastomer1	349	1.74	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_45	\N	easy_hard_axis	Ni2MnGa_sample_1	45	1	\N	0
strain_sample1_elastomer1_tens_1_11	\N	strain	sample1_elastomer1	11	0.05	\N	0
strain_wire_actuator_NiTi#6_110	\N	strain	wire_actuator_NiTi#6	110	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_135	\N	mechanic_stress	sample1_elastomer1	135	493300	\N	0
strain_sample1_elastomer1_tens_1_31	\N	strain	sample1_elastomer1	31	0.15	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_47	\N	easy_hard_axis	Ni2MnGa_sample_1	47	1	\N	0
strain_sample1_elastomer1_tens_1_136	\N	strain	sample1_elastomer1	136	0.675	\N	0
mechanic_stress_sample1_elastomer1_tens_1_306	\N	mechanic_stress	sample1_elastomer1	306	975300	\N	0
temperature_wire_actuator_1_133	\N	temperature	wire_actuator_1	133	50	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_126	\N	magnetic_field_strength	Ni2MnGa_sample_1	126	250000	\N	0
strain_sample1_elastomer1_tens_1_295	\N	strain	sample1_elastomer1	295	1.47	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_44	\N	easy_hard_axis	Ni2MnGa_sample_1	44	1	\N	0
strain_sample1_elastomer1_tens_1_293	\N	strain	sample1_elastomer1	293	1.46	\N	0
hysteresis_wire_actuator_1_135	\N	hysteresis	wire_actuator_1	135	1	\N	0
strain_sample1_elastomer1_tens_1_163	\N	strain	sample1_elastomer1	163	0.81	\N	0
temperature_wire_actuator_1_120	\N	temperature	wire_actuator_1	120	71	\N	0
M420_1_10	\N	C00	M420	0	1.4166399999999998e-08	\N	0
\N	\N	vacuum_permeability	Ni2MnGa_sample_1	0	1.2566e-06	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_13	\N	magnetic_field_strength	Ni2MnGa_sample_1	13	120000	\N	0
strain_sample1_elastomer1_tens_1_95	\N	strain	sample1_elastomer1	95	0.47	\N	0
strain_sample1_elastomer1_tens_1_23	\N	strain	sample1_elastomer1	23	0.11	\N	0
strain_sample1_elastomer1_tens_1_174	\N	strain	sample1_elastomer1	174	0.865	\N	0
strain_sample1_elastomer1_tens_1_38	\N	strain	sample1_elastomer1	38	0.185	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_140	\N	easy_hard_axis	Ni2MnGa_sample_1	140	0	\N	0
M420_1_7	\N	B00	M420	0	3.55e-10	\N	0
strain_sample1_elastomer1_tens_1_299	\N	strain	sample1_elastomer1	299	1.49	\N	0
strain_sample1_elastomer1_tens_1_159	\N	strain	sample1_elastomer1	159	0.79	\N	0
mechanic_stress_sample1_elastomer1_tens_1_321	\N	mechanic_stress	sample1_elastomer1	321	1028300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_246	\N	mechanic_stress	sample1_elastomer1	246	789300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_17	\N	magnetic_field_strength	Ni2MnGa_sample_1	17	160000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_224	\N	mechanic_stress	sample1_elastomer1	224	726300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_123	\N	mechanic_stress	sample1_elastomer1	123	463300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_149	\N	easy_hard_axis	Ni2MnGa_sample_1	149	0	\N	0
strain_sample1_elastomer1_tens_1_149	\N	strain	sample1_elastomer1	149	0.74	\N	0
strain_sample1_elastomer1_tens_1_344	\N	strain	sample1_elastomer1	344	1.715	\N	0
strain_wire_actuator_NiTi#6_219	\N	strain	wire_actuator_NiTi#6	219	0.045	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_24	\N	magnetic_field_strength	Ni2MnGa_sample_1	24	230000	\N	0
hysteresis_wire_actuator_1_137	\N	hysteresis	wire_actuator_1	137	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_260	\N	mechanic_stress	sample1_elastomer1	260	830300	\N	0
strain_wire_actuator_NiTi#6_214	\N	strain	wire_actuator_NiTi#6	214	0.02	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_122	\N	magnetic_field_strength	Ni2MnGa_sample_1	122	210000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_264	\N	mechanic_stress	sample1_elastomer1	264	842300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_22	\N	mechanic_stress	sample1_elastomer1	22	118300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_142	\N	mechanic_stress	sample1_elastomer1	142	512300	\N	0
strain_sample1_elastomer1_tens_1_389	\N	strain	sample1_elastomer1	389	1.94	\N	0
mechanic_stress_sample1_elastomer1_tens_1_195	\N	mechanic_stress	sample1_elastomer1	195	649300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_38	\N	mechanic_stress	sample1_elastomer1	38	191300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_275	\N	mechanic_stress	sample1_elastomer1	275	876300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_39	\N	mechanic_stress	sample1_elastomer1	39	196300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_13	\N	mechanic_stress	sample1_elastomer1	13	71000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_110	\N	magnetic_field_strength	Ni2MnGa_sample_1	110	90000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_120	\N	easy_hard_axis	Ni2MnGa_sample_1	120	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_68	\N	mechanic_stress	sample1_elastomer1	68	303300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_108	\N	mechanic_stress	sample1_elastomer1	108	422300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_44	\N	mechanic_stress	sample1_elastomer1	44	216300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_144	\N	easy_hard_axis	Ni2MnGa_sample_1	144	0	\N	0
strain_sample1_elastomer1_tens_1_50	\N	strain	sample1_elastomer1	50	0.245	\N	0
mechanic_stress_sample1_elastomer1_tens_1_61	\N	mechanic_stress	sample1_elastomer1	61	281300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_280	\N	mechanic_stress	sample1_elastomer1	280	891800	\N	0
strain_sample1_elastomer1_tens_1_20	\N	strain	sample1_elastomer1	20	0.095	\N	0
mechanic_stress_sample1_elastomer1_tens_1_350	\N	mechanic_stress	sample1_elastomer1	350	1123300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_117	\N	mechanic_stress	sample1_elastomer1	117	446300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_235	\N	mechanic_stress	sample1_elastomer1	235	756300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_373	\N	mechanic_stress	sample1_elastomer1	373	1198300	\N	0
hysteresis_wire_actuator_1_121	\N	hysteresis	wire_actuator_1	121	2	\N	0
mechanic_stress_sample1_elastomer1_tens_1_375	\N	mechanic_stress	sample1_elastomer1	375	1208300	\N	0
strain_sample1_elastomer1_tens_1_379	\N	strain	sample1_elastomer1	379	1.89	\N	0
mechanic_stress_sample1_elastomer1_tens_1_401	\N	mechanic_stress	sample1_elastomer1	401	1308300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_362	\N	mechanic_stress	sample1_elastomer1	362	1163300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_259	\N	mechanic_stress	sample1_elastomer1	259	827300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_3	\N	easy_hard_axis	Ni2MnGa_sample_1	3	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_57	\N	magnetic_field_strength	Ni2MnGa_sample_1	57	560000	\N	0
strain_sample1_elastomer1_tens_1_220	\N	strain	sample1_elastomer1	220	1.095	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_108	\N	easy_hard_axis	Ni2MnGa_sample_1	108	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_361	\N	mechanic_stress	sample1_elastomer1	361	1158300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_146	\N	easy_hard_axis	Ni2MnGa_sample_1	146	0	\N	0
strain_sample1_elastomer1_tens_1_262	\N	strain	sample1_elastomer1	262	1.305	\N	0
mechanic_stress_sample1_elastomer1_tens_1_239	\N	mechanic_stress	sample1_elastomer1	239	768300	\N	0
strain_sample1_elastomer1_tens_1_4	\N	strain	sample1_elastomer1	4	0.015	\N	0
strain_sample1_elastomer1_tens_1_265	\N	strain	sample1_elastomer1	265	1.32	\N	0
strain_sample1_elastomer1_tens_1_237	\N	strain	sample1_elastomer1	237	1.18	\N	0
mechanic_stress_sample1_elastomer1_tens_1_168	\N	mechanic_stress	sample1_elastomer1	168	578300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_380	\N	mechanic_stress	sample1_elastomer1	380	1223300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_26	\N	easy_hard_axis	Ni2MnGa_sample_1	26	1	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_131	\N	easy_hard_axis	Ni2MnGa_sample_1	131	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_400	\N	mechanic_stress	sample1_elastomer1	400	1303300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_115	\N	mechanic_stress	sample1_elastomer1	115	441300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_117	\N	magnetic_field_strength	Ni2MnGa_sample_1	117	160000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_42	\N	easy_hard_axis	Ni2MnGa_sample_1	42	1	\N	0
martensite_content_wire_actuator_NiTi#6_318	\N	martensite_content	wire_actuator_NiTi#6	318	1	\N	0
twinning_stress_stick_Ni2MnGa_sample_1	\N	twinning_stress	stick_Ni2MnGa_sample_1	0	300000	\N	0
strain_sample1_elastomer1_tens_1_277	\N	strain	sample1_elastomer1	277	1.38	\N	0
mechanic_stress_sample1_elastomer1_tens_1_258	\N	mechanic_stress	sample1_elastomer1	258	824300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_207	\N	mechanic_stress	sample1_elastomer1	207	679300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_113	\N	mechanic_stress	sample1_elastomer1	113	435300	\N	0
martensite_content_wire_actuator_NiTi#6_314	\N	martensite_content	wire_actuator_NiTi#6	314	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_101	\N	magnetic_field_strength	Ni2MnGa_sample_1	101	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_166	\N	mechanic_stress	sample1_elastomer1	166	572300	\N	0
strain_sample1_elastomer1_tens_1_348	\N	strain	sample1_elastomer1	348	1.735	\N	0
strain_sample1_elastomer1_tens_1_85	\N	strain	sample1_elastomer1	85	0.42	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_29	\N	easy_hard_axis	Ni2MnGa_sample_1	29	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_147	\N	magnetic_field_strength	Ni2MnGa_sample_1	147	460000	\N	0
temperature_wire_actuator_1_112	\N	temperature	wire_actuator_1	112	48	\N	0
strain_sample1_elastomer1_tens_1_224	\N	strain	sample1_elastomer1	224	1.115	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_111	\N	easy_hard_axis	Ni2MnGa_sample_1	111	0	\N	0
strain_sample1_elastomer1_tens_1_100	\N	strain	sample1_elastomer1	100	0.495	\N	0
mechanic_stress_sample1_elastomer1_tens_1_86	\N	mechanic_stress	sample1_elastomer1	86	359300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_150	\N	magnetic_field_strength	Ni2MnGa_sample_1	150	490000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_345	\N	mechanic_stress	sample1_elastomer1	345	1108300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_34	\N	easy_hard_axis	Ni2MnGa_sample_1	34	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_323	\N	mechanic_stress	sample1_elastomer1	323	1028300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_379	\N	mechanic_stress	sample1_elastomer1	379	1218300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_161	\N	mechanic_stress	sample1_elastomer1	161	560300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_98	\N	mechanic_stress	sample1_elastomer1	98	394300	\N	0
martensite_finish_temperature_wire_actuator_1_0	\N	martensite_finish_temperature	wire_actuator_1	0	47.2	\N	0
mechanic_stress_sample1_elastomer1_tens_1_181	\N	mechanic_stress	sample1_elastomer1	181	612300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_74	\N	mechanic_stress	sample1_elastomer1	74	323300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_53	\N	mechanic_stress	sample1_elastomer1	53	250300	\N	0
\N	\N	vacuum_permeability	Ni2MnGa	0	1.2566e-06	\N	0
mechanic_stress_sample1_elastomer1_tens_1_152	\N	mechanic_stress	sample1_elastomer1	152	538300	\N	0
strain_sample1_elastomer1_tens_1_385	\N	strain	sample1_elastomer1	385	1.92	\N	0
strain_sample1_elastomer1_tens_1_91	\N	strain	sample1_elastomer1	91	0.45	\N	0
mechanic_stress_sample1_elastomer1_tens_1_149	\N	mechanic_stress	sample1_elastomer1	149	530300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_16	\N	magnetic_field_strength	Ni2MnGa_sample_1	16	150000	\N	0
strain_sample1_elastomer1_tens_1_19	\N	strain	sample1_elastomer1	19	0.09	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_125	\N	magnetic_field_strength	Ni2MnGa_sample_1	125	240000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_31	\N	easy_hard_axis	Ni2MnGa_sample_1	31	1	\N	0
strain_sample1_elastomer1_tens_1_357	\N	strain	sample1_elastomer1	357	1.78	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_126	\N	easy_hard_axis	Ni2MnGa_sample_1	126	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_211	\N	mechanic_stress	sample1_elastomer1	211	691300	\N	0
strain_sample1_elastomer1_tens_1_16	\N	strain	sample1_elastomer1	16	0.075	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_5	\N	easy_hard_axis	Ni2MnGa_sample_1	5	1	\N	0
strain_sample1_elastomer1_tens_1_51	\N	strain	sample1_elastomer1	51	0.25	\N	0
temperature_wire_actuator_1_141	\N	temperature	wire_actuator_1	141	72	\N	0
strain_sample1_elastomer1_tens_1_304	\N	strain	sample1_elastomer1	304	1.515	\N	0
mechanic_stress_sample1_elastomer1_tens_1_218	\N	mechanic_stress	sample1_elastomer1	218	709800	\N	0
strain_sample1_elastomer1_tens_1_121	\N	strain	sample1_elastomer1	121	0.6	\N	0
mechanic_stress_sample1_elastomer1_tens_1_111	\N	mechanic_stress	sample1_elastomer1	111	430300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_157	\N	magnetic_field_strength	Ni2MnGa_sample_1	157	560000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_15	\N	magnetic_field_strength	Ni2MnGa_sample_1	15	140000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_217	\N	mechanic_stress	sample1_elastomer1	217	707300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_41	\N	magnetic_field_strength	Ni2MnGa_sample_1	41	400000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_31	\N	mechanic_stress	sample1_elastomer1	31	161300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_387	\N	mechanic_stress	sample1_elastomer1	387	1248300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_231	\N	mechanic_stress	sample1_elastomer1	231	745300	\N	0
strain_sample1_elastomer1_tens_1_56	\N	strain	sample1_elastomer1	56	0.275	\N	0
strain_wire_actuator_NiTi#6_213	\N	strain	wire_actuator_NiTi#6	213	0.015	\N	0
mechanic_stress_sample1_elastomer1_tens_1_9	\N	mechanic_stress	sample1_elastomer1	9	47333	\N	0
mechanic_stress_sample1_elastomer1_tens_1_285	\N	mechanic_stress	sample1_elastomer1	285	907300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_8	\N	magnetic_field_strength	Ni2MnGa_sample_1	8	70000	\N	0
strain_sample1_elastomer1_tens_1_329	\N	strain	sample1_elastomer1	329	1.64	\N	0
mechanic_stress_sample1_elastomer1_tens_1_62	\N	mechanic_stress	sample1_elastomer1	62	283300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_359	\N	mechanic_stress	sample1_elastomer1	359	1158300	\N	0
M420_1_4	\N	A00	M420	0	-5.7e-12	\N	0
mechanic_stress_sample1_elastomer1_tens_1_206	\N	mechanic_stress	sample1_elastomer1	206	676300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_12	\N	magnetic_field_strength	Ni2MnGa_sample_1	12	110000	\N	0
martensite_content_wire_actuator_NiTi#6_118	\N	martensite_content	wire_actuator_NiTi#6	118	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_177	\N	mechanic_stress	sample1_elastomer1	177	602300	\N	0
strain_wire_actuator_NiTi#6_222	\N	strain	wire_actuator_NiTi#6	222	0.06	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_147	\N	easy_hard_axis	Ni2MnGa_sample_1	147	0	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_133	\N	easy_hard_axis	Ni2MnGa_sample_1	133	0	\N	0
strain_sample1_elastomer1_tens_1_383	\N	strain	sample1_elastomer1	383	1.91	\N	0
mechanic_stress_sample1_elastomer1_tens_1_105	\N	mechanic_stress	sample1_elastomer1	105	414300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_210	\N	mechanic_stress	sample1_elastomer1	210	688300	\N	0
strain_sample1_elastomer1_tens_1_40	\N	strain	sample1_elastomer1	40	0.195	\N	0
mechanic_stress_sample1_elastomer1_tens_1_140	\N	mechanic_stress	sample1_elastomer1	140	506300	\N	0
strain_sample1_elastomer1_tens_1_394	\N	strain	sample1_elastomer1	394	1.965	\N	0
strain_sample1_elastomer1_tens_1_284	\N	strain	sample1_elastomer1	284	1.415	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_102	\N	easy_hard_axis	Ni2MnGa_sample_1	102	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_312	\N	mechanic_stress	sample1_elastomer1	312	993300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_125	\N	mechanic_stress	sample1_elastomer1	125	468300	\N	0
hysteresis_wire_actuator_1_117	\N	hysteresis	wire_actuator_1	117	2	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_148	\N	magnetic_field_strength	Ni2MnGa_sample_1	148	470000	\N	0
strain_wire_actuator_NiTi#6_220	\N	strain	wire_actuator_NiTi#6	220	0.05	\N	0
strain_wire_actuator_NiTi#6_316	\N	strain	wire_actuator_NiTi#6	316	0.03	\N	0
\N	\N	vacuum_permittivity	NiTi	0	8.854e-12	\N	0
strain_sample1_elastomer1_tens_1_156	\N	strain	sample1_elastomer1	156	0.775	\N	0
strain_sample1_elastomer1_tens_1_202	\N	strain	sample1_elastomer1	202	1.005	\N	0
mechanic_stress_sample1_elastomer1_tens_1_366	\N	mechanic_stress	sample1_elastomer1	366	1178300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_129	\N	mechanic_stress	sample1_elastomer1	129	477300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_137	\N	easy_hard_axis	Ni2MnGa_sample_1	137	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_51	\N	mechanic_stress	sample1_elastomer1	51	243300	\N	0
strain_sample1_elastomer1_tens_1_103	\N	strain	sample1_elastomer1	103	0.51	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_106	\N	magnetic_field_strength	Ni2MnGa_sample_1	106	50000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_26	\N	magnetic_field_strength	Ni2MnGa_sample_1	26	250000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_347	\N	mechanic_stress	sample1_elastomer1	347	1118300	\N	0
hysteresis_wire_actuator_1_120	\N	hysteresis	wire_actuator_1	120	2	\N	0
temperature_wire_actuator_1_136	\N	temperature	wire_actuator_1	136	53	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_158	\N	easy_hard_axis	Ni2MnGa_sample_1	158	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_172	\N	mechanic_stress	sample1_elastomer1	172	589300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_377	\N	mechanic_stress	sample1_elastomer1	377	1218300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_122	\N	mechanic_stress	sample1_elastomer1	122	460300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_118	\N	mechanic_stress	sample1_elastomer1	118	448300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_381	\N	mechanic_stress	sample1_elastomer1	381	1228300	\N	0
hysteresis_wire_actuator_1_119	\N	hysteresis	wire_actuator_1	119	2	\N	0
strain_sample1_elastomer1_tens_1_185	\N	strain	sample1_elastomer1	185	0.92	\N	0
martensite_content_wire_actuator_NiTi#6_320	\N	martensite_content	wire_actuator_NiTi#6	320	1	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_101	\N	easy_hard_axis	Ni2MnGa_sample_1	101	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_8	\N	mechanic_stress	sample1_elastomer1	8	41533	\N	0
mechanic_stress_sample1_elastomer1_tens_1_299	\N	mechanic_stress	sample1_elastomer1	299	952300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_191	\N	mechanic_stress	sample1_elastomer1	191	638300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_334	\N	mechanic_stress	sample1_elastomer1	334	1073300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_14	\N	mechanic_stress	sample1_elastomer1	14	76033	\N	0
strain_sample1_elastomer1_tens_1_322	\N	strain	sample1_elastomer1	322	1.605	\N	0
mechanic_stress_sample1_elastomer1_tens_1_89	\N	mechanic_stress	sample1_elastomer1	89	369300	\N	0
strain_sample1_elastomer1_tens_1_388	\N	strain	sample1_elastomer1	388	1.935	\N	0
strain_sample1_elastomer1_tens_1_157	\N	strain	sample1_elastomer1	157	0.78	\N	0
strain_sample1_elastomer1_tens_1_271	\N	strain	sample1_elastomer1	271	1.35	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_132	\N	easy_hard_axis	Ni2MnGa_sample_1	132	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_241	\N	mechanic_stress	sample1_elastomer1	241	774300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_4	\N	magnetic_field_strength	Ni2MnGa_sample_1	4	30000	\N	0
strain_sample1_elastomer1_tens_1_218	\N	strain	sample1_elastomer1	218	1.085	\N	0
mechanic_stress_sample1_elastomer1_tens_1_352	\N	mechanic_stress	sample1_elastomer1	352	1128300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_249	\N	mechanic_stress	sample1_elastomer1	249	798300	\N	0
M420_1_9	\N	C00	M420	0	1.4166399999999998e-08	\N	0
mechanic_stress_sample1_elastomer1_tens_1_283	\N	mechanic_stress	sample1_elastomer1	283	901300	\N	0
\N	\N	vacuum_permittivity	Ni2MnGa	0	8.854e-12	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_47	\N	magnetic_field_strength	Ni2MnGa_sample_1	47	460000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_230	\N	mechanic_stress	sample1_elastomer1	230	742800	\N	0
strain_sample1_elastomer1_tens_1_243	\N	strain	sample1_elastomer1	243	1.21	\N	0
mechanic_stress_sample1_elastomer1_tens_1_160	\N	mechanic_stress	sample1_elastomer1	160	558300	\N	0
strain_sample1_elastomer1_tens_1_36	\N	strain	sample1_elastomer1	36	0.175	\N	0
strain_sample1_elastomer1_tens_1_361	\N	strain	sample1_elastomer1	361	1.8	\N	0
strain_sample1_elastomer1_tens_1_61	\N	strain	sample1_elastomer1	61	0.3	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_135	\N	magnetic_field_strength	Ni2MnGa_sample_1	135	340000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_176	\N	mechanic_stress	sample1_elastomer1	176	598300	\N	0
strain_sample1_elastomer1_tens_1_235	\N	strain	sample1_elastomer1	235	1.17	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_105	\N	easy_hard_axis	Ni2MnGa_sample_1	105	0	\N	0
strain_sample1_elastomer1_tens_1_144	\N	strain	sample1_elastomer1	144	0.715	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_160	\N	easy_hard_axis	Ni2MnGa_sample_1	160	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_26	\N	mechanic_stress	sample1_elastomer1	26	138300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_162	\N	easy_hard_axis	Ni2MnGa_sample_1	162	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_226	\N	mechanic_stress	sample1_elastomer1	226	732300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_128	\N	easy_hard_axis	Ni2MnGa_sample_1	128	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_17	\N	mechanic_stress	sample1_elastomer1	17	92300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_6	\N	mechanic_stress	sample1_elastomer1	6	30400	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_107	\N	magnetic_field_strength	Ni2MnGa_sample_1	107	60000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_92	\N	mechanic_stress	sample1_elastomer1	92	378300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_40	\N	easy_hard_axis	Ni2MnGa_sample_1	40	1	\N	0
strain_sample1_elastomer1_tens_1_197	\N	strain	sample1_elastomer1	197	0.98	\N	0
strain_sample1_elastomer1_tens_1_5	\N	strain	sample1_elastomer1	5	0.02	\N	0
mechanic_stress_sample1_elastomer1_tens_1_236	\N	mechanic_stress	sample1_elastomer1	236	759800	\N	0
strain_sample1_elastomer1_tens_1_122	\N	strain	sample1_elastomer1	122	0.605	\N	0
strain_sample1_elastomer1_tens_1_105	\N	strain	sample1_elastomer1	105	0.52	\N	0
strain_sample1_elastomer1_tens_1_3	\N	strain	sample1_elastomer1	3	0.01	\N	0
strain_sample1_elastomer1_tens_1_192	\N	strain	sample1_elastomer1	192	0.955	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_159	\N	magnetic_field_strength	Ni2MnGa_sample_1	159	580000	\N	0
temperature_wire_actuator_1_142	\N	temperature	wire_actuator_1	142	73.5	\N	0
strain_sample1_elastomer1_tens_1_362	\N	strain	sample1_elastomer1	362	1.805	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_16	\N	easy_hard_axis	Ni2MnGa_sample_1	16	1	\N	0
strain_sample1_elastomer1_tens_1_347	\N	strain	sample1_elastomer1	347	1.73	\N	0
mechanic_stress_sample1_elastomer1_tens_1_294	\N	mechanic_stress	sample1_elastomer1	294	936300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_106	\N	mechanic_stress	sample1_elastomer1	106	416300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_23	\N	magnetic_field_strength	Ni2MnGa_sample_1	23	220000	\N	0
strain_sample1_elastomer1_tens_1_39	\N	strain	sample1_elastomer1	39	0.19	\N	0
mechanic_stress_sample1_elastomer1_tens_1_76	\N	mechanic_stress	sample1_elastomer1	76	329300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_57	\N	easy_hard_axis	Ni2MnGa_sample_1	57	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_180	\N	mechanic_stress	sample1_elastomer1	180	610300	\N	0
strain_sample1_elastomer1_tens_1_326	\N	strain	sample1_elastomer1	326	1.625	\N	0
max_block_stress_M420	\N	maximal_blocking_stress	M420	0	1907689.37	\N	0
mechanic_stress_sample1_elastomer1_tens_1_7	\N	mechanic_stress	sample1_elastomer1	7	35900	\N	0
mechanic_stress_sample1_elastomer1_tens_1_242	\N	mechanic_stress	sample1_elastomer1	242	776300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_47	\N	mechanic_stress	sample1_elastomer1	47	229300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_119	\N	easy_hard_axis	Ni2MnGa_sample_1	119	0	\N	0
strain_sample1_elastomer1_tens_1_363	\N	strain	sample1_elastomer1	363	1.81	\N	0
mechanic_stress_sample1_elastomer1_tens_1_300	\N	mechanic_stress	sample1_elastomer1	300	955800	\N	0
temperature_wire_actuator_1_119	\N	temperature	wire_actuator_1	119	70	\N	0
mechanic_stress_sample1_elastomer1_tens_1_356	\N	mechanic_stress	sample1_elastomer1	356	1143300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_55	\N	mechanic_stress	sample1_elastomer1	55	257300	\N	0
strain_sample1_elastomer1_tens_1_356	\N	strain	sample1_elastomer1	356	1.775	\N	0
mechanic_stress_sample1_elastomer1_tens_1_10	\N	mechanic_stress	sample1_elastomer1	10	52900	\N	0
mechanic_stress_sample1_elastomer1_tens_1_386	\N	mechanic_stress	sample1_elastomer1	386	1248300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_337	\N	mechanic_stress	sample1_elastomer1	337	1078300	\N	0
strain_sample1_elastomer1_tens_1_71	\N	strain	sample1_elastomer1	71	0.35	\N	0
strain_sample1_elastomer1_tens_1_135	\N	strain	sample1_elastomer1	135	0.67	\N	0
hysteresis_wire_actuator_1_122	\N	hysteresis	wire_actuator_1	122	2	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_128	\N	magnetic_field_strength	Ni2MnGa_sample_1	128	270000	\N	0
strain_sample1_elastomer1_tens_1_195	\N	strain	sample1_elastomer1	195	0.97	\N	0
martensite_content_wire_actuator_NiTi#6_114	\N	martensite_content	wire_actuator_NiTi#6	114	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_331	\N	mechanic_stress	sample1_elastomer1	331	1058300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_3	\N	magnetic_field_strength	Ni2MnGa_sample_1	3	20000	\N	0
strain_sample1_elastomer1_tens_1_172	\N	strain	sample1_elastomer1	172	0.855	\N	0
strain_sample1_elastomer1_tens_1_34	\N	strain	sample1_elastomer1	34	0.165	\N	0
mechanic_stress_sample1_elastomer1_tens_1_344	\N	mechanic_stress	sample1_elastomer1	344	1103300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_355	\N	mechanic_stress	sample1_elastomer1	355	1138300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_49	\N	magnetic_field_strength	Ni2MnGa_sample_1	49	480000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_338	\N	mechanic_stress	sample1_elastomer1	338	1083300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_4	\N	mechanic_stress	sample1_elastomer1	4	17200	\N	0
martensite_content_wire_actuator_NiTi#6_311	\N	martensite_content	wire_actuator_NiTi#6	311	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_281	\N	mechanic_stress	sample1_elastomer1	281	894300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_169	\N	mechanic_stress	sample1_elastomer1	169	581300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_147	\N	mechanic_stress	sample1_elastomer1	147	525300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_97	\N	mechanic_stress	sample1_elastomer1	97	391300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_358	\N	mechanic_stress	sample1_elastomer1	358	1153300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_201	\N	mechanic_stress	sample1_elastomer1	201	665300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_40	\N	magnetic_field_strength	Ni2MnGa_sample_1	40	390000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_134	\N	easy_hard_axis	Ni2MnGa_sample_1	134	0	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_23	\N	easy_hard_axis	Ni2MnGa_sample_1	23	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_53	\N	magnetic_field_strength	Ni2MnGa_sample_1	53	520000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_378	\N	mechanic_stress	sample1_elastomer1	378	1218300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_35	\N	easy_hard_axis	Ni2MnGa_sample_1	35	1	\N	0
strain_sample1_elastomer1_tens_1_212	\N	strain	sample1_elastomer1	212	1.055	\N	0
martensite_content_wire_actuator_NiTi#6_119	\N	martensite_content	wire_actuator_NiTi#6	119	0	\N	0
strain_sample1_elastomer1_tens_1_114	\N	strain	sample1_elastomer1	114	0.565	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_140	\N	magnetic_field_strength	Ni2MnGa_sample_1	140	390000	\N	0
strain_sample1_elastomer1_tens_1_116	\N	strain	sample1_elastomer1	116	0.575	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_115	\N	magnetic_field_strength	Ni2MnGa_sample_1	115	140000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_49	\N	mechanic_stress	sample1_elastomer1	49	235300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_392	\N	mechanic_stress	sample1_elastomer1	392	1268300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_371	\N	mechanic_stress	sample1_elastomer1	371	1198300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_149	\N	magnetic_field_strength	Ni2MnGa_sample_1	149	480000	\N	0
strain_sample1_elastomer1_tens_1_213	\N	strain	sample1_elastomer1	213	1.06	\N	0
strain_sample1_elastomer1_tens_1_253	\N	strain	sample1_elastomer1	253	1.26	\N	0
mechanic_stress_sample1_elastomer1_tens_1_291	\N	mechanic_stress	sample1_elastomer1	291	927300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_198	\N	mechanic_stress	sample1_elastomer1	198	656300	\N	0
M420_2_1	\N	B01	M420	0	-0.011	\N	0
strain_sample1_elastomer1_tens_1_398	\N	strain	sample1_elastomer1	398	1.985	\N	0
mechanic_stress_sample1_elastomer1_tens_1_340	\N	mechanic_stress	sample1_elastomer1	340	1088300	\N	0
\N	\N	vacuum_permittivity	Ni2MnGa_sample_1	0	8.854e-12	\N	0
mechanic_stress_sample1_elastomer1_tens_1_360	\N	mechanic_stress	sample1_elastomer1	360	1158300	\N	0
martensite_content_wire_actuator_NiTi#6_322	\N	martensite_content	wire_actuator_NiTi#6	322	1	\N	0
relative_permittivity_el_1	\N	relative_permittivity	elastomer_1	0	3.2	\N	0
mechanic_stress_sample1_elastomer1_tens_1_390	\N	mechanic_stress	sample1_elastomer1	390	1263300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_393	\N	mechanic_stress	sample1_elastomer1	393	1268300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_121	\N	easy_hard_axis	Ni2MnGa_sample_1	121	0	\N	0
strain_sample1_elastomer1_tens_1_252	\N	strain	sample1_elastomer1	252	1.255	\N	0
strain_sample1_elastomer1_tens_1_134	\N	strain	sample1_elastomer1	134	0.665	\N	0
strain_sample1_elastomer1_tens_1_193	\N	strain	sample1_elastomer1	193	0.96	\N	0
mechanic_stress_sample1_elastomer1_tens_1_368	\N	mechanic_stress	sample1_elastomer1	368	1183300	\N	0
strain_sample1_elastomer1_tens_1_77	\N	strain	sample1_elastomer1	77	0.38	\N	0
mechanic_stress_sample1_elastomer1_tens_1_252	\N	mechanic_stress	sample1_elastomer1	252	806800	\N	0
strain_sample1_elastomer1_tens_1_247	\N	strain	sample1_elastomer1	247	1.23	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_142	\N	easy_hard_axis	Ni2MnGa_sample_1	142	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_170	\N	mechanic_stress	sample1_elastomer1	170	584300	\N	0
strain_sample1_elastomer1_tens_1_67	\N	strain	sample1_elastomer1	67	0.33	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_152	\N	easy_hard_axis	Ni2MnGa_sample_1	152	0	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_154	\N	magnetic_field_strength	Ni2MnGa_sample_1	154	530000	\N	0
flux_density_Ni2MnGa_sample1_meas_flux_over_field_2_101	\N	flux_density	Ni2MnGa_sample_1	101	0	\N	0
strain_sample1_elastomer1_tens_1_314	\N	strain	sample1_elastomer1	314	1.565	\N	0
mechanic_stress_sample1_elastomer1_tens_1_27	\N	mechanic_stress	sample1_elastomer1	27	143300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_173	\N	mechanic_stress	sample1_elastomer1	173	591300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_77	\N	mechanic_stress	sample1_elastomer1	77	333300	\N	0
strain_sample1_elastomer1_tens_1_32	\N	strain	sample1_elastomer1	32	0.155	\N	0
martensite_content_wire_actuator_NiTi#6_222	\N	martensite_content	wire_actuator_NiTi#6	222	0.5	\N	0
mechanic_stress_sample1_elastomer1_tens_1_215	\N	mechanic_stress	sample1_elastomer1	215	701300	\N	0
strain_sample1_elastomer1_tens_1_345	\N	strain	sample1_elastomer1	345	1.72	\N	0
strain_sample1_elastomer1_tens_1_98	\N	strain	sample1_elastomer1	98	0.485	\N	0
strain_sample1_elastomer1_tens_1_89	\N	strain	sample1_elastomer1	89	0.44	\N	0
mechanic_stress_sample1_elastomer1_tens_1_131	\N	mechanic_stress	sample1_elastomer1	131	483300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_34	\N	mechanic_stress	sample1_elastomer1	34	174300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_84	\N	mechanic_stress	sample1_elastomer1	84	353300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_48	\N	easy_hard_axis	Ni2MnGa_sample_1	48	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_25	\N	mechanic_stress	sample1_elastomer1	25	133300	\N	0
M420_1_8	\N	B00	M420	0	5.25e-10	\N	0
martensite_content_wire_actuator_NiTi#6_321	\N	martensite_content	wire_actuator_NiTi#6	321	1	\N	0
strain_sample1_elastomer1_tens_1_400	\N	strain	sample1_elastomer1	400	1.995	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_58	\N	easy_hard_axis	Ni2MnGa_sample_1	58	1	\N	0
\N	\N	vacuum_permeability	sample1_elastomer1	0	1.2566e-06	\N	0
strain_sample1_elastomer1_tens_1_338	\N	strain	sample1_elastomer1	338	1.685	\N	0
mechanic_stress_sample1_elastomer1_tens_1_398	\N	mechanic_stress	sample1_elastomer1	398	1293300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_9	\N	magnetic_field_strength	Ni2MnGa_sample_1	9	80000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_223	\N	mechanic_stress	sample1_elastomer1	223	723300	\N	0
strain_sample1_elastomer1_tens_1_328	\N	strain	sample1_elastomer1	328	1.635	\N	0
strain_sample1_elastomer1_tens_1_2	\N	strain	sample1_elastomer1	2	0.005	\N	0
strain_sample1_elastomer1_tens_1_170	\N	strain	sample1_elastomer1	170	0.845	\N	0
mechanic_stress_sample1_elastomer1_tens_1_388	\N	mechanic_stress	sample1_elastomer1	388	1253300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_27	\N	magnetic_field_strength	Ni2MnGa_sample_1	27	260000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_33	\N	easy_hard_axis	Ni2MnGa_sample_1	33	1	\N	0
strain_sample1_elastomer1_tens_1_275	\N	strain	sample1_elastomer1	275	1.37	\N	0
strain_sample1_elastomer1_tens_1_324	\N	strain	sample1_elastomer1	324	1.615	\N	0
strain_sample1_elastomer1_tens_1_315	\N	strain	sample1_elastomer1	315	1.57	\N	0
strain_sample1_elastomer1_tens_1_303	\N	strain	sample1_elastomer1	303	1.51	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_144	\N	magnetic_field_strength	Ni2MnGa_sample_1	144	430000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_50	\N	magnetic_field_strength	Ni2MnGa_sample_1	50	490000	\N	0
temperature_wire_actuator_1_117	\N	temperature	wire_actuator_1	117	60	\N	0
mechanic_stress_sample1_elastomer1_tens_1_102	\N	mechanic_stress	sample1_elastomer1	102	406300	\N	0
strain_sample1_elastomer1_tens_1_132	\N	strain	sample1_elastomer1	132	0.655	\N	0
mechanic_stress_sample1_elastomer1_tens_1_104	\N	mechanic_stress	sample1_elastomer1	104	410300	\N	0
strain_sample1_elastomer1_tens_1_52	\N	strain	sample1_elastomer1	52	0.255	\N	0
strain_sample1_elastomer1_tens_1_268	\N	strain	sample1_elastomer1	268	1.335	\N	0
martensite_content_wire_actuator_NiTi#6_211	\N	martensite_content	wire_actuator_NiTi#6	211	0.5	\N	0
strain_sample1_elastomer1_tens_1_334	\N	strain	sample1_elastomer1	334	1.665	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_46	\N	easy_hard_axis	Ni2MnGa_sample_1	46	1	\N	0
strain_sample1_elastomer1_tens_1_207	\N	strain	sample1_elastomer1	207	1.03	\N	0
mechanic_stress_sample1_elastomer1_tens_1_229	\N	mechanic_stress	sample1_elastomer1	229	740300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_82	\N	mechanic_stress	sample1_elastomer1	82	348300	\N	0
hysteresis_wire_actuator_1_114	\N	hysteresis	wire_actuator_1	114	2	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_124	\N	easy_hard_axis	Ni2MnGa_sample_1	124	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_328	\N	mechanic_stress	sample1_elastomer1	328	1048300	\N	0
strain_sample1_elastomer1_tens_1_111	\N	strain	sample1_elastomer1	111	0.55	\N	0
strain_sample1_elastomer1_tens_1_321	\N	strain	sample1_elastomer1	321	1.6	\N	0
strain_sample1_elastomer1_tens_1_101	\N	strain	sample1_elastomer1	101	0.5	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_52	\N	magnetic_field_strength	Ni2MnGa_sample_1	52	510000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_41	\N	easy_hard_axis	Ni2MnGa_sample_1	41	1	\N	0
dielectric_strength_el_1	\N	dielectric_strength	elastomer_1	0	99000000	\N	0
strain_sample1_elastomer1_tens_1_8	\N	strain	sample1_elastomer1	8	0.035	\N	0
strain_wire_actuator_NiTi#6_115	\N	strain	wire_actuator_NiTi#6	115	0.025	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_49	\N	easy_hard_axis	Ni2MnGa_sample_1	49	1	\N	0
strain_sample1_elastomer1_tens_1_263	\N	strain	sample1_elastomer1	263	1.31	\N	0
hysteresis_wire_actuator_1_112	\N	hysteresis	wire_actuator_1	112	2	\N	0
martensite_content_wire_actuator_NiTi#6_220	\N	martensite_content	wire_actuator_NiTi#6	220	0.5	\N	0
mechanic_stress_sample1_elastomer1_tens_1_137	\N	mechanic_stress	sample1_elastomer1	137	499300	\N	0
strain_wire_actuator_NiTi#6_312	\N	strain	wire_actuator_NiTi#6	312	0.01	\N	0
strain_sample1_elastomer1_tens_1_261	\N	strain	sample1_elastomer1	261	1.3	\N	0
mechanic_stress_sample1_elastomer1_tens_1_56	\N	mechanic_stress	sample1_elastomer1	56	262300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_269	\N	mechanic_stress	sample1_elastomer1	269	857300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_2	\N	mechanic_stress	sample1_elastomer1	2	5600	\N	0
mechanic_stress_sample1_elastomer1_tens_1_227	\N	mechanic_stress	sample1_elastomer1	227	735300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_268	\N	mechanic_stress	sample1_elastomer1	268	854300	\N	0
temperature_wire_actuator_1_113	\N	temperature	wire_actuator_1	113	50	\N	0
strain_sample1_elastomer1_tens_1_341	\N	strain	sample1_elastomer1	341	1.7	\N	0
mechanic_stress_sample1_elastomer1_tens_1_254	\N	mechanic_stress	sample1_elastomer1	254	812800	\N	0
mechanic_stress_sample1_elastomer1_tens_1_24	\N	mechanic_stress	sample1_elastomer1	24	129300	\N	0
strain_sample1_elastomer1_tens_1_387	\N	strain	sample1_elastomer1	387	1.93	\N	0
strain_sample1_elastomer1_tens_1_166	\N	strain	sample1_elastomer1	166	0.825	\N	0
mechanic_stress_sample1_elastomer1_tens_1_233	\N	mechanic_stress	sample1_elastomer1	233	750300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_16	\N	mechanic_stress	sample1_elastomer1	16	87300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_145	\N	magnetic_field_strength	Ni2MnGa_sample_1	145	440000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_6	\N	easy_hard_axis	Ni2MnGa_sample_1	6	1	\N	0
strain_sample1_elastomer1_tens_1_294	\N	strain	sample1_elastomer1	294	1.465	\N	0
mechanic_stress_sample1_elastomer1_tens_1_57	\N	mechanic_stress	sample1_elastomer1	57	265300	\N	0
strain_sample1_elastomer1_tens_1_201	\N	strain	sample1_elastomer1	201	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_232	\N	mechanic_stress	sample1_elastomer1	232	747800	\N	0
martensite_content_wire_actuator_NiTi#6_213	\N	martensite_content	wire_actuator_NiTi#6	213	0.5	\N	0
strain_sample1_elastomer1_tens_1_18	\N	strain	sample1_elastomer1	18	0.085	\N	0
mechanic_stress_sample1_elastomer1_tens_1_96	\N	mechanic_stress	sample1_elastomer1	96	388300	\N	0
flux_density_Ni2MnGa_sample1_meas_flux_over_field_1_1	\N	flux_density	Ni2MnGa_sample_1	1	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_148	\N	mechanic_stress	sample1_elastomer1	148	527300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_148	\N	easy_hard_axis	Ni2MnGa_sample_1	148	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_42	\N	mechanic_stress	sample1_elastomer1	42	208300	\N	0
strain_wire_actuator_NiTi#6_321	\N	strain	wire_actuator_NiTi#6	321	0.055	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_59	\N	magnetic_field_strength	Ni2MnGa_sample_1	59	580000	\N	0
strain_sample1_elastomer1_tens_1_272	\N	strain	sample1_elastomer1	272	1.355	\N	0
mechanic_stress_sample1_elastomer1_tens_1_103	\N	mechanic_stress	sample1_elastomer1	103	409300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_24	\N	easy_hard_axis	Ni2MnGa_sample_1	24	1	\N	0
M420_1_1	\N	A00	M420	0	1.54e-11	\N	0
mechanic_stress_sample1_elastomer1_tens_1_200	\N	mechanic_stress	sample1_elastomer1	200	663300	\N	0
strain_sample1_elastomer1_tens_1_279	\N	strain	sample1_elastomer1	279	1.39	\N	0
strain_sample1_elastomer1_tens_1_10	\N	strain	sample1_elastomer1	10	0.045	\N	0
austenite_finish_temperature_wire_actuator_1_0	\N	austenite_finish_temperature	wire_actuator_1	0	73.9	\N	0
strain_sample1_elastomer1_tens_1_205	\N	strain	sample1_elastomer1	205	1.02	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_160	\N	magnetic_field_strength	Ni2MnGa_sample_1	160	590000	\N	0
strain_wire_actuator_NiTi#6_113	\N	strain	wire_actuator_NiTi#6	113	0.015	\N	0
mechanic_stress_sample1_elastomer1_tens_1_385	\N	mechanic_stress	sample1_elastomer1	385	1248300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_127	\N	magnetic_field_strength	Ni2MnGa_sample_1	127	260000	\N	0
strain_sample1_elastomer1_tens_1_242	\N	strain	sample1_elastomer1	242	1.205	\N	0
strain_sample1_elastomer1_tens_1_300	\N	strain	sample1_elastomer1	300	1.495	\N	0
hysteresis_wire_actuator_1_142	\N	hysteresis	wire_actuator_1	142	1	\N	0
strain_wire_actuator_NiTi#6_217	\N	strain	wire_actuator_NiTi#6	217	0.035	\N	0
strain_sample1_elastomer1_tens_1_179	\N	strain	sample1_elastomer1	179	0.89	\N	0
mechanic_stress_sample1_elastomer1_tens_1_134	\N	mechanic_stress	sample1_elastomer1	134	491300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_341	\N	mechanic_stress	sample1_elastomer1	341	1088300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_150	\N	easy_hard_axis	Ni2MnGa_sample_1	150	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_155	\N	mechanic_stress	sample1_elastomer1	155	545300	\N	0
strain_sample1_elastomer1_tens_1_392	\N	strain	sample1_elastomer1	392	1.955	\N	0
martensite_content_wire_actuator_NiTi#6_112	\N	martensite_content	wire_actuator_NiTi#6	112	0	\N	0
strain_sample1_elastomer1_tens_1_119	\N	strain	sample1_elastomer1	119	0.59	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_43	\N	magnetic_field_strength	Ni2MnGa_sample_1	43	420000	\N	0
strain_sample1_elastomer1_tens_1_42	\N	strain	sample1_elastomer1	42	0.205	\N	0
strain_sample1_elastomer1_tens_1_66	\N	strain	sample1_elastomer1	66	0.325	\N	0
mechanic_stress_sample1_elastomer1_tens_1_186	\N	mechanic_stress	sample1_elastomer1	186	625300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_72	\N	mechanic_stress	sample1_elastomer1	72	316300	\N	0
strain_sample1_elastomer1_tens_1_270	\N	strain	sample1_elastomer1	270	1.345	\N	0
mechanic_stress_sample1_elastomer1_tens_1_95	\N	mechanic_stress	sample1_elastomer1	95	385300	\N	0
strain_wire_actuator_NiTi#6_119	\N	strain	wire_actuator_NiTi#6	119	0.045	\N	0
mechanic_stress_sample1_elastomer1_tens_1_19	\N	mechanic_stress	sample1_elastomer1	19	103300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_42	\N	magnetic_field_strength	Ni2MnGa_sample_1	42	410000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_274	\N	mechanic_stress	sample1_elastomer1	274	873800	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_6	\N	magnetic_field_strength	Ni2MnGa_sample_1	6	50000	\N	0
strain_sample1_elastomer1_tens_1_211	\N	strain	sample1_elastomer1	211	1.05	\N	0
strain_sample1_elastomer1_tens_1_373	\N	strain	sample1_elastomer1	373	1.86	\N	0
mechanic_stress_sample1_elastomer1_tens_1_343	\N	mechanic_stress	sample1_elastomer1	343	1098300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_313	\N	mechanic_stress	sample1_elastomer1	313	998300	\N	0
strain_sample1_elastomer1_tens_1_158	\N	strain	sample1_elastomer1	158	0.785	\N	0
strain_sample1_elastomer1_tens_1_140	\N	strain	sample1_elastomer1	140	0.695	\N	0
strain_sample1_elastomer1_tens_1_257	\N	strain	sample1_elastomer1	257	1.28	\N	0
strain_sample1_elastomer1_tens_1_276	\N	strain	sample1_elastomer1	276	1.375	\N	0
mechanic_stress_sample1_elastomer1_tens_1_309	\N	mechanic_stress	sample1_elastomer1	309	988300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_62	\N	easy_hard_axis	Ni2MnGa_sample_1	62	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_190	\N	mechanic_stress	sample1_elastomer1	190	636300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_189	\N	mechanic_stress	sample1_elastomer1	189	633300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_270	\N	mechanic_stress	sample1_elastomer1	270	860800	\N	0
strain_sample1_elastomer1_tens_1_186	\N	strain	sample1_elastomer1	186	0.925	\N	0
strain_sample1_elastomer1_tens_1_310	\N	strain	sample1_elastomer1	310	1.545	\N	0
mechanic_stress_sample1_elastomer1_tens_1_369	\N	mechanic_stress	sample1_elastomer1	369	1188300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_316	\N	mechanic_stress	sample1_elastomer1	316	1008300	\N	0
temperature_wire_actuator_1_115	\N	temperature	wire_actuator_1	115	52	\N	0
temperature_wire_actuator_1_143	\N	temperature	wire_actuator_1	143	75	\N	0
strain_sample1_elastomer1_tens_1_305	\N	strain	sample1_elastomer1	305	1.52	\N	0
strain_sample1_elastomer1_tens_1_274	\N	strain	sample1_elastomer1	274	1.365	\N	0
mechanic_stress_sample1_elastomer1_tens_1_276	\N	mechanic_stress	sample1_elastomer1	276	878800	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_141	\N	magnetic_field_strength	Ni2MnGa_sample_1	141	400000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_63	\N	mechanic_stress	sample1_elastomer1	63	287300	\N	0
strain_sample1_elastomer1_tens_1_229	\N	strain	sample1_elastomer1	229	1.14	\N	0
plateau_start_strain_wire_actuator_NiTi#6_0	\N	plateau_start_strain	wire_actuator_NiTi#6	0	0.0088	\N	0
mechanic_stress_sample1_elastomer1_tens_1_157	\N	mechanic_stress	sample1_elastomer1	157	550300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_138	\N	magnetic_field_strength	Ni2MnGa_sample_1	138	370000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_228	\N	mechanic_stress	sample1_elastomer1	228	737800	\N	0
temperature_wire_actuator_1_116	\N	temperature	wire_actuator_1	116	53	\N	0
mechanic_stress_sample1_elastomer1_tens_1_32	\N	mechanic_stress	sample1_elastomer1	32	166300	\N	0
M420_2_2	\N	B01	M420	0	0.025	\N	0
mechanic_stress_sample1_elastomer1_tens_1_272	\N	mechanic_stress	sample1_elastomer1	272	867800	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_153	\N	easy_hard_axis	Ni2MnGa_sample_1	153	0	\N	0
lattice_constant_c_Ni2MnGa_sample_1	\N	lattice_constant_c	Ni2MnGa_sample_1	0	0.55	\N	0
strain_sample1_elastomer1_tens_1_396	\N	strain	sample1_elastomer1	396	1.975	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_157	\N	easy_hard_axis	Ni2MnGa_sample_1	157	0	\N	0
strain_sample1_elastomer1_tens_1_371	\N	strain	sample1_elastomer1	371	1.85	\N	0
strain_sample1_elastomer1_tens_1_162	\N	strain	sample1_elastomer1	162	0.805	\N	0
strain_sample1_elastomer1_tens_1_171	\N	strain	sample1_elastomer1	171	0.85	\N	0
strain_sample1_elastomer1_tens_1_297	\N	strain	sample1_elastomer1	297	1.48	\N	0
strain_sample1_elastomer1_tens_1_291	\N	strain	sample1_elastomer1	291	1.45	\N	0
mechanic_stress_sample1_elastomer1_tens_1_144	\N	mechanic_stress	sample1_elastomer1	144	517300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_297	\N	mechanic_stress	sample1_elastomer1	297	945300	\N	0
strain_sample1_elastomer1_tens_1_313	\N	strain	sample1_elastomer1	313	1.56	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_135	\N	easy_hard_axis	Ni2MnGa_sample_1	135	0	\N	0
temperature_wire_actuator_1_132	\N	temperature	wire_actuator_1	132	48	\N	0
mechanic_stress_sample1_elastomer1_tens_1_179	\N	mechanic_stress	sample1_elastomer1	179	607300	\N	0
strain_sample1_elastomer1_tens_1_245	\N	strain	sample1_elastomer1	245	1.22	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_136	\N	easy_hard_axis	Ni2MnGa_sample_1	136	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_159	\N	mechanic_stress	sample1_elastomer1	159	555300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_146	\N	mechanic_stress	sample1_elastomer1	146	522300	\N	0
strain_sample1_elastomer1_tens_1_368	\N	strain	sample1_elastomer1	368	1.835	\N	0
strain_sample1_elastomer1_tens_1_110	\N	strain	sample1_elastomer1	110	0.545	\N	0
temperature_wire_actuator_1_138	\N	temperature	wire_actuator_1	138	65	\N	0
mechanic_stress_sample1_elastomer1_tens_1_310	\N	mechanic_stress	sample1_elastomer1	310	988300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_41	\N	mechanic_stress	sample1_elastomer1	41	204300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_54	\N	magnetic_field_strength	Ni2MnGa_sample_1	54	530000	\N	0
strain_sample1_elastomer1_tens_1_60	\N	strain	sample1_elastomer1	60	0.295	\N	0
strain_sample1_elastomer1_tens_1_264	\N	strain	sample1_elastomer1	264	1.315	\N	0
mechanic_stress_sample1_elastomer1_tens_1_112	\N	mechanic_stress	sample1_elastomer1	112	433300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_15	\N	mechanic_stress	sample1_elastomer1	15	82300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_5	\N	magnetic_field_strength	Ni2MnGa_sample_1	5	40000	\N	0
M420_1_2	\N	A00	M420	0	1.87e-11	\N	0
strain_sample1_elastomer1_tens_1_308	\N	strain	sample1_elastomer1	308	1.535	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_19	\N	magnetic_field_strength	Ni2MnGa_sample_1	19	180000	\N	0
strain_sample1_elastomer1_tens_1_81	\N	strain	sample1_elastomer1	81	0.4	\N	0
mechanic_stress_sample1_elastomer1_tens_1_307	\N	mechanic_stress	sample1_elastomer1	307	978300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_165	\N	mechanic_stress	sample1_elastomer1	165	570300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_152	\N	magnetic_field_strength	Ni2MnGa_sample_1	152	510000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_25	\N	magnetic_field_strength	Ni2MnGa_sample_1	25	240000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_143	\N	easy_hard_axis	Ni2MnGa_sample_1	143	0	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_27	\N	easy_hard_axis	Ni2MnGa_sample_1	27	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_56	\N	magnetic_field_strength	Ni2MnGa_sample_1	56	550000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_278	\N	mechanic_stress	sample1_elastomer1	278	885300	\N	0
strain_sample1_elastomer1_tens_1_358	\N	strain	sample1_elastomer1	358	1.785	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_38	\N	magnetic_field_strength	Ni2MnGa_sample_1	38	370000	\N	0
strain_wire_actuator_NiTi#6_320	\N	strain	wire_actuator_NiTi#6	320	0.05	\N	0
mechanic_stress_sample1_elastomer1_tens_1_213	\N	mechanic_stress	sample1_elastomer1	213	696300	\N	0
strain_sample1_elastomer1_tens_1_167	\N	strain	sample1_elastomer1	167	0.83	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_112	\N	magnetic_field_strength	Ni2MnGa_sample_1	112	110000	\N	0
strain_wire_actuator_NiTi#6_211	\N	strain	wire_actuator_NiTi#6	211	0.005	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_33	\N	magnetic_field_strength	Ni2MnGa_sample_1	33	320000	\N	0
strain_sample1_elastomer1_tens_1_120	\N	strain	sample1_elastomer1	120	0.595	\N	0
mechanic_stress_sample1_elastomer1_tens_1_48	\N	mechanic_stress	sample1_elastomer1	48	231300	\N	0
strain_sample1_elastomer1_tens_1_176	\N	strain	sample1_elastomer1	176	0.875	\N	0
mechanic_stress_sample1_elastomer1_tens_1_204	\N	mechanic_stress	sample1_elastomer1	204	670800	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_20	\N	easy_hard_axis	Ni2MnGa_sample_1	20	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_193	\N	mechanic_stress	sample1_elastomer1	193	643300	\N	0
strain_sample1_elastomer1_tens_1_59	\N	strain	sample1_elastomer1	59	0.29	\N	0
strain_sample1_elastomer1_tens_1_188	\N	strain	sample1_elastomer1	188	0.935	\N	0
strain_sample1_elastomer1_tens_1_301	\N	strain	sample1_elastomer1	301	1.5	\N	0
strain_sample1_elastomer1_tens_1_258	\N	strain	sample1_elastomer1	258	1.285	\N	0
strain_wire_actuator_NiTi#6_319	\N	strain	wire_actuator_NiTi#6	319	0.045	\N	0
hysteresis_wire_actuator_1_113	\N	hysteresis	wire_actuator_1	113	2	\N	0
strain_sample1_elastomer1_tens_1_350	\N	strain	sample1_elastomer1	350	1.745	\N	0
mechanic_stress_sample1_elastomer1_tens_1_364	\N	mechanic_stress	sample1_elastomer1	364	1173300	\N	0
strain_sample1_elastomer1_tens_1_164	\N	strain	sample1_elastomer1	164	0.815	\N	0
mechanic_stress_sample1_elastomer1_tens_1_67	\N	mechanic_stress	sample1_elastomer1	67	300300	\N	0
strain_sample1_elastomer1_tens_1_364	\N	strain	sample1_elastomer1	364	1.815	\N	0
mechanic_stress_sample1_elastomer1_tens_1_326	\N	mechanic_stress	sample1_elastomer1	326	1043300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_99	\N	mechanic_stress	sample1_elastomer1	99	397300	\N	0
martensite_content_wire_actuator_NiTi#6_215	\N	martensite_content	wire_actuator_NiTi#6	215	0.5	\N	0
strain_sample1_elastomer1_tens_1_287	\N	strain	sample1_elastomer1	287	1.43	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_129	\N	magnetic_field_strength	Ni2MnGa_sample_1	129	280000	\N	0
M420_1_5	\N	A00	M420	0	4.5e-11	\N	0
strain_sample1_elastomer1_tens_1_191	\N	strain	sample1_elastomer1	191	0.95	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_22	\N	easy_hard_axis	Ni2MnGa_sample_1	22	1	\N	0
temperature_wire_actuator_1_134	\N	temperature	wire_actuator_1	134	51	\N	0
strain_sample1_elastomer1_tens_1_160	\N	strain	sample1_elastomer1	160	0.795	\N	0
strain_sample1_elastomer1_tens_1_204	\N	strain	sample1_elastomer1	204	1.015	\N	0
strain_sample1_elastomer1_tens_1_65	\N	strain	sample1_elastomer1	65	0.32	\N	0
hysteresis_wire_actuator_1_141	\N	hysteresis	wire_actuator_1	141	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_262	\N	mechanic_stress	sample1_elastomer1	262	836300	\N	0
initial_strain_elastomer_1	\N	initial_strain	sample1_elastomer1	0	0.05	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_22	\N	magnetic_field_strength	Ni2MnGa_sample_1	22	210000	\N	0
martensite_content_wire_actuator_NiTi#6_219	\N	martensite_content	wire_actuator_NiTi#6	219	0.5	\N	0
strain_sample1_elastomer1_tens_1_376	\N	strain	sample1_elastomer1	376	1.875	\N	0
strain_sample1_elastomer1_tens_1_79	\N	strain	sample1_elastomer1	79	0.39	\N	0
strain_sample1_elastomer1_tens_1_26	\N	strain	sample1_elastomer1	26	0.125	\N	0
mechanic_stress_sample1_elastomer1_tens_1_208	\N	mechanic_stress	sample1_elastomer1	208	682300	\N	0
martensite_content_wire_actuator_NiTi#6_212	\N	martensite_content	wire_actuator_NiTi#6	212	0.5	\N	0
strain_sample1_elastomer1_tens_1_127	\N	strain	sample1_elastomer1	127	0.63	\N	0
strain_sample1_elastomer1_tens_1_106	\N	strain	sample1_elastomer1	106	0.525	\N	0
martensite_content_wire_actuator_NiTi#6_111	\N	martensite_content	wire_actuator_NiTi#6	111	0	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_114	\N	magnetic_field_strength	Ni2MnGa_sample_1	114	130000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_71	\N	mechanic_stress	sample1_elastomer1	71	313300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_105	\N	magnetic_field_strength	Ni2MnGa_sample_1	105	40000	\N	0
martensite_content_wire_actuator_NiTi#6_115	\N	martensite_content	wire_actuator_NiTi#6	115	0	\N	0
martensite_content_wire_actuator_NiTi#6_500	\N	martensite_content	wire_actuator_NiTi#6	500	1	\N	0
strain_sample1_elastomer1_tens_1_255	\N	strain	sample1_elastomer1	255	1.27	\N	0
mechanic_stress_sample1_elastomer1_tens_1_205	\N	mechanic_stress	sample1_elastomer1	205	673300	\N	0
strain_sample1_elastomer1_tens_1_346	\N	strain	sample1_elastomer1	346	1.725	\N	0
strain_sample1_elastomer1_tens_1_246	\N	strain	sample1_elastomer1	246	1.225	\N	0
mechanic_stress_sample1_elastomer1_tens_1_153	\N	mechanic_stress	sample1_elastomer1	153	539300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_336	\N	mechanic_stress	sample1_elastomer1	336	1078300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_125	\N	easy_hard_axis	Ni2MnGa_sample_1	125	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_302	\N	mechanic_stress	sample1_elastomer1	302	962300	\N	0
strain_sample1_elastomer1_tens_1_316	\N	strain	sample1_elastomer1	316	1.575	\N	0
strain_wire_actuator_NiTi#6_210	\N	strain	wire_actuator_NiTi#6	210	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_332	\N	mechanic_stress	sample1_elastomer1	332	1063300	\N	0
strain_sample1_elastomer1_tens_1_292	\N	strain	sample1_elastomer1	292	1.455	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_106	\N	easy_hard_axis	Ni2MnGa_sample_1	106	0	\N	0
strain_sample1_elastomer1_tens_1_254	\N	strain	sample1_elastomer1	254	1.265	\N	0
martensite_content_wire_actuator_NiTi#6_312	\N	martensite_content	wire_actuator_NiTi#6	312	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_123	\N	magnetic_field_strength	Ni2MnGa_sample_1	123	220000	\N	0
strain_sample1_elastomer1_tens_1_128	\N	strain	sample1_elastomer1	128	0.635	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_50	\N	easy_hard_axis	Ni2MnGa_sample_1	50	1	\N	0
strain_sample1_elastomer1_tens_1_17	\N	strain	sample1_elastomer1	17	0.08	\N	0
strain_sample1_elastomer1_tens_1_27	\N	strain	sample1_elastomer1	27	0.13	\N	0
mechanic_stress_sample1_elastomer1_tens_1_290	\N	mechanic_stress	sample1_elastomer1	290	923800	\N	0
strain_sample1_elastomer1_tens_1_72	\N	strain	sample1_elastomer1	72	0.355	\N	0
strain_wire_actuator_NiTi#6_117	\N	strain	wire_actuator_NiTi#6	117	0.035	\N	0
mechanic_stress_sample1_elastomer1_tens_1_322	\N	mechanic_stress	sample1_elastomer1	322	1028300	\N	0
strain_sample1_elastomer1_tens_1_189	\N	strain	sample1_elastomer1	189	0.94	\N	0
strain_sample1_elastomer1_tens_1_93	\N	strain	sample1_elastomer1	93	0.46	\N	0
mechanic_stress_sample1_elastomer1_tens_1_1	\N	mechanic_stress	sample1_elastomer1	1	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_363	\N	mechanic_stress	sample1_elastomer1	363	1168300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_151	\N	mechanic_stress	sample1_elastomer1	151	535300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_107	\N	mechanic_stress	sample1_elastomer1	107	419300	\N	0
lattice_constant_a_Ni2MnGa_sample_1	\N	lattice_constant_a	Ni2MnGa_sample_1	0	0.59	\N	0
strain_sample1_elastomer1_tens_1_343	\N	strain	sample1_elastomer1	343	1.71	\N	0
mechanic_stress_sample1_elastomer1_tens_1_265	\N	mechanic_stress	sample1_elastomer1	265	845300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_194	\N	mechanic_stress	sample1_elastomer1	194	646300	\N	0
strain_sample1_elastomer1_tens_1_141	\N	strain	sample1_elastomer1	141	0.7	\N	0
mechanic_stress_sample1_elastomer1_tens_1_109	\N	mechanic_stress	sample1_elastomer1	109	423300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_295	\N	mechanic_stress	sample1_elastomer1	295	939300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_120	\N	magnetic_field_strength	Ni2MnGa_sample_1	120	190000	\N	0
strain_sample1_elastomer1_tens_1_401	\N	strain	sample1_elastomer1	401	2	\N	0
strain_sample1_elastomer1_tens_1_251	\N	strain	sample1_elastomer1	251	1.25	\N	0
strain_sample1_elastomer1_tens_1_151	\N	strain	sample1_elastomer1	151	0.75	\N	0
strain_sample1_elastomer1_tens_1_280	\N	strain	sample1_elastomer1	280	1.395	\N	0
strain_sample1_elastomer1_tens_1_206	\N	strain	sample1_elastomer1	206	1.025	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_31	\N	magnetic_field_strength	Ni2MnGa_sample_1	31	300000	\N	0
strain_sample1_elastomer1_tens_1_281	\N	strain	sample1_elastomer1	281	1.4	\N	0
mechanic_stress_sample1_elastomer1_tens_1_351	\N	mechanic_stress	sample1_elastomer1	351	1128300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_346	\N	mechanic_stress	sample1_elastomer1	346	1113300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_15	\N	easy_hard_axis	Ni2MnGa_sample_1	15	1	\N	0
strain_sample1_elastomer1_tens_1_6	\N	strain	sample1_elastomer1	6	0.025	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_58	\N	magnetic_field_strength	Ni2MnGa_sample_1	58	570000	\N	0
E-modulus_twinned_martensite_wire_actuator_NiTi#6_0	\N	E-modulus_twinned_martensite	wire_actuator_NiTi#6	0	22800000000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_35	\N	mechanic_stress	sample1_elastomer1	35	178300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_21	\N	easy_hard_axis	Ni2MnGa_sample_1	21	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_102	\N	magnetic_field_strength	Ni2MnGa_sample_1	102	10000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_136	\N	magnetic_field_strength	Ni2MnGa_sample_1	136	350000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_119	\N	mechanic_stress	sample1_elastomer1	119	452300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_126	\N	mechanic_stress	sample1_elastomer1	126	471300	\N	0
strain_wire_actuator_NiTi#6_120	\N	strain	wire_actuator_NiTi#6	120	0.05	\N	0
strain_sample1_elastomer1_tens_1_131	\N	strain	sample1_elastomer1	131	0.65	\N	0
strain_sample1_elastomer1_tens_1_35	\N	strain	sample1_elastomer1	35	0.17	\N	0
mechanic_stress_sample1_elastomer1_tens_1_327	\N	mechanic_stress	sample1_elastomer1	327	1048300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_374	\N	mechanic_stress	sample1_elastomer1	374	1203300	\N	0
strain_sample1_elastomer1_tens_1_282	\N	strain	sample1_elastomer1	282	1.405	\N	0
strain_sample1_elastomer1_tens_1_340	\N	strain	sample1_elastomer1	340	1.695	\N	0
\N	\N	vacuum_permeability	elastomer_1	0	1.2566e-06	\N	0
strain_sample1_elastomer1_tens_1_53	\N	strain	sample1_elastomer1	53	0.26	\N	0
hysteresis_wire_actuator_1_140	\N	hysteresis	wire_actuator_1	140	1	\N	0
strain_sample1_elastomer1_tens_1_90	\N	strain	sample1_elastomer1	90	0.445	\N	0
strain_sample1_elastomer1_tens_1_104	\N	strain	sample1_elastomer1	104	0.515	\N	0
strain_sample1_elastomer1_tens_1_80	\N	strain	sample1_elastomer1	80	0.395	\N	0
strain_sample1_elastomer1_tens_1_259	\N	strain	sample1_elastomer1	259	1.29	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_30	\N	magnetic_field_strength	Ni2MnGa_sample_1	30	290000	\N	0
strain_sample1_elastomer1_tens_1_75	\N	strain	sample1_elastomer1	75	0.37	\N	0
mechanic_stress_sample1_elastomer1_tens_1_255	\N	mechanic_stress	sample1_elastomer1	255	815300	\N	0
temperature_wire_actuator_1_122	\N	temperature	wire_actuator_1	122	73.5	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_110	\N	easy_hard_axis	Ni2MnGa_sample_1	110	0	\N	0
strain_wire_actuator_NiTi#6_121	\N	strain	wire_actuator_NiTi#6	121	0.055	\N	0
mechanic_stress_sample1_elastomer1_tens_1_83	\N	mechanic_stress	sample1_elastomer1	83	351300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_70	\N	mechanic_stress	sample1_elastomer1	70	310300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_93	\N	mechanic_stress	sample1_elastomer1	93	380300	\N	0
strain_wire_actuator_NiTi#6_111	\N	strain	wire_actuator_NiTi#6	111	0.005	\N	0
strain_sample1_elastomer1_tens_1_336	\N	strain	sample1_elastomer1	336	1.675	\N	0
strain_sample1_elastomer1_tens_1_133	\N	strain	sample1_elastomer1	133	0.66	\N	0
mechanic_stress_sample1_elastomer1_tens_1_234	\N	mechanic_stress	sample1_elastomer1	234	753300	\N	0
strain_wire_actuator_NiTi#6_310	\N	strain	wire_actuator_NiTi#6	310	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_143	\N	mechanic_stress	sample1_elastomer1	143	513300	\N	0
strain_sample1_elastomer1_tens_1_47	\N	strain	sample1_elastomer1	47	0.23	\N	0
hardening_parameter_wire_actuator_NiTi#6_0	\N	hardening_parameter	wire_actuator_NiTi#6	0	150000000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_39	\N	easy_hard_axis	Ni2MnGa_sample_1	39	1	\N	0
strain_sample1_elastomer1_tens_1_146	\N	strain	sample1_elastomer1	146	0.725	\N	0
strain_sample1_elastomer1_tens_1_86	\N	strain	sample1_elastomer1	86	0.425	\N	0
strain_sample1_elastomer1_tens_1_73	\N	strain	sample1_elastomer1	73	0.36	\N	0
strain_sample1_elastomer1_tens_1_88	\N	strain	sample1_elastomer1	88	0.435	\N	0
strain_sample1_elastomer1_tens_1_76	\N	strain	sample1_elastomer1	76	0.375	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_55	\N	easy_hard_axis	Ni2MnGa_sample_1	55	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_261	\N	mechanic_stress	sample1_elastomer1	261	833300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_127	\N	easy_hard_axis	Ni2MnGa_sample_1	127	0	\N	0
strain_sample1_elastomer1_tens_1_331	\N	strain	sample1_elastomer1	331	1.65	\N	0
temperature_wire_actuator_1_118	\N	temperature	wire_actuator_1	118	65	\N	0
temperature_wire_actuator_1_121	\N	temperature	wire_actuator_1	121	72	\N	0
mechanic_stress_sample1_elastomer1_tens_1_257	\N	mechanic_stress	sample1_elastomer1	257	821300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_73	\N	mechanic_stress	sample1_elastomer1	73	319300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_25	\N	easy_hard_axis	Ni2MnGa_sample_1	25	1	\N	0
strain_sample1_elastomer1_tens_1_94	\N	strain	sample1_elastomer1	94	0.465	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_116	\N	magnetic_field_strength	Ni2MnGa_sample_1	116	150000	\N	0
strain_sample1_elastomer1_tens_1_175	\N	strain	sample1_elastomer1	175	0.87	\N	0
strain_sample1_elastomer1_tens_1_14	\N	strain	sample1_elastomer1	14	0.065	\N	0
strain_sample1_elastomer1_tens_1_30	\N	strain	sample1_elastomer1	30	0.145	\N	0
mechanic_stress_sample1_elastomer1_tens_1_277	\N	mechanic_stress	sample1_elastomer1	277	881300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_131	\N	magnetic_field_strength	Ni2MnGa_sample_1	131	300000	\N	0
strain_sample1_elastomer1_tens_1_69	\N	strain	sample1_elastomer1	69	0.34	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_37	\N	magnetic_field_strength	Ni2MnGa_sample_1	37	360000	\N	0
strain_sample1_elastomer1_tens_1_397	\N	strain	sample1_elastomer1	397	1.98	\N	0
mechanic_stress_sample1_elastomer1_tens_1_271	\N	mechanic_stress	sample1_elastomer1	271	864300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_248	\N	mechanic_stress	sample1_elastomer1	248	795300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_104	\N	magnetic_field_strength	Ni2MnGa_sample_1	104	30000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_244	\N	mechanic_stress	sample1_elastomer1	244	782300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_94	\N	mechanic_stress	sample1_elastomer1	94	383300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_8	\N	easy_hard_axis	Ni2MnGa_sample_1	8	1	\N	0
strain_sample1_elastomer1_tens_1_366	\N	strain	sample1_elastomer1	366	1.825	\N	0
strain_sample1_elastomer1_tens_1_113	\N	strain	sample1_elastomer1	113	0.56	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_159	\N	easy_hard_axis	Ni2MnGa_sample_1	159	0	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_19	\N	easy_hard_axis	Ni2MnGa_sample_1	19	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_175	\N	mechanic_stress	sample1_elastomer1	175	596300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_5	\N	mechanic_stress	sample1_elastomer1	5	24567	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_53	\N	easy_hard_axis	Ni2MnGa_sample_1	53	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_39	\N	magnetic_field_strength	Ni2MnGa_sample_1	39	380000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_151	\N	magnetic_field_strength	Ni2MnGa_sample_1	151	500000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_113	\N	magnetic_field_strength	Ni2MnGa_sample_1	113	120000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_90	\N	mechanic_stress	sample1_elastomer1	90	371300	\N	0
strain_sample1_elastomer1_tens_1_369	\N	strain	sample1_elastomer1	369	1.84	\N	0
strain_sample1_elastomer1_tens_1_25	\N	strain	sample1_elastomer1	25	0.12	\N	0
strain_sample1_elastomer1_tens_1_267	\N	strain	sample1_elastomer1	267	1.33	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_13	\N	easy_hard_axis	Ni2MnGa_sample_1	13	1	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_113	\N	easy_hard_axis	Ni2MnGa_sample_1	113	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_397	\N	mechanic_stress	sample1_elastomer1	397	1288300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_145	\N	mechanic_stress	sample1_elastomer1	145	519300	\N	0
strain_sample1_elastomer1_tens_1_273	\N	strain	sample1_elastomer1	273	1.36	\N	0
mechanic_stress_sample1_elastomer1_tens_1_305	\N	mechanic_stress	sample1_elastomer1	305	972300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_116	\N	easy_hard_axis	Ni2MnGa_sample_1	116	0	\N	0
strain_sample1_elastomer1_tens_1_102	\N	strain	sample1_elastomer1	102	0.505	\N	0
mechanic_stress_sample1_elastomer1_tens_1_45	\N	mechanic_stress	sample1_elastomer1	45	220300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_107	\N	easy_hard_axis	Ni2MnGa_sample_1	107	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_150	\N	mechanic_stress	sample1_elastomer1	150	532300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_365	\N	mechanic_stress	sample1_elastomer1	365	1178300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_196	\N	mechanic_stress	sample1_elastomer1	196	651300	\N	0
martensite_content_wire_actuator_NiTi#6_122	\N	martensite_content	wire_actuator_NiTi#6	122	0	\N	0
hysteresis_wire_actuator_1_115	\N	hysteresis	wire_actuator_1	115	2	\N	0
strain_sample1_elastomer1_tens_1_382	\N	strain	sample1_elastomer1	382	1.905	\N	0
strain_sample1_elastomer1_tens_1_325	\N	strain	sample1_elastomer1	325	1.62	\N	0
strain_sample1_elastomer1_tens_1_317	\N	strain	sample1_elastomer1	317	1.58	\N	0
strain_sample1_elastomer1_tens_1_260	\N	strain	sample1_elastomer1	260	1.295	\N	0
strain_sample1_elastomer1_tens_1_194	\N	strain	sample1_elastomer1	194	0.965	\N	0
hysteresis_wire_actuator_1_139	\N	hysteresis	wire_actuator_1	139	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_282	\N	mechanic_stress	sample1_elastomer1	282	897800	\N	0
mechanic_stress_sample1_elastomer1_tens_1_212	\N	mechanic_stress	sample1_elastomer1	212	693800	\N	0
strain_sample1_elastomer1_tens_1_109	\N	strain	sample1_elastomer1	109	0.54	\N	0
mechanic_stress_sample1_elastomer1_tens_1_81	\N	mechanic_stress	sample1_elastomer1	81	344300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_118	\N	easy_hard_axis	Ni2MnGa_sample_1	118	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_40	\N	mechanic_stress	sample1_elastomer1	40	200300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_156	\N	magnetic_field_strength	Ni2MnGa_sample_1	156	550000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_367	\N	mechanic_stress	sample1_elastomer1	367	1178300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_253	\N	mechanic_stress	sample1_elastomer1	253	810300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_163	\N	mechanic_stress	sample1_elastomer1	163	565300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_51	\N	magnetic_field_strength	Ni2MnGa_sample_1	51	500000	\N	0
strain_sample1_elastomer1_tens_1_234	\N	strain	sample1_elastomer1	234	1.165	\N	0
strain_sample1_elastomer1_tens_1_391	\N	strain	sample1_elastomer1	391	1.95	\N	0
mechanic_stress_sample1_elastomer1_tens_1_263	\N	mechanic_stress	sample1_elastomer1	263	839300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_132	\N	mechanic_stress	sample1_elastomer1	132	485300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_1	\N	magnetic_field_strength	Ni2MnGa_sample_1	1	0	\N	0
strain_sample1_elastomer1_tens_1_68	\N	strain	sample1_elastomer1	68	0.335	\N	0
strain_sample1_elastomer1_tens_1_214	\N	strain	sample1_elastomer1	214	1.065	\N	0
strain_sample1_elastomer1_tens_1_168	\N	strain	sample1_elastomer1	168	0.835	\N	0
strain_sample1_elastomer1_tens_1_57	\N	strain	sample1_elastomer1	57	0.28	\N	0
strain_sample1_elastomer1_tens_1_152	\N	strain	sample1_elastomer1	152	0.755	\N	0
strain_sample1_elastomer1_tens_1_296	\N	strain	sample1_elastomer1	296	1.475	\N	0
mechanic_stress_sample1_elastomer1_tens_1_287	\N	mechanic_stress	sample1_elastomer1	287	913300	\N	0
strain_sample1_elastomer1_tens_1_318	\N	strain	sample1_elastomer1	318	1.585	\N	0
mechanic_stress_sample1_elastomer1_tens_1_138	\N	mechanic_stress	sample1_elastomer1	138	502300	\N	0
hysteresis_wire_actuator_1_111	\N	hysteresis	wire_actuator_1	111	2	\N	0
mechanic_stress_sample1_elastomer1_tens_1_203	\N	mechanic_stress	sample1_elastomer1	203	668300	\N	0
strain_sample1_elastomer1_tens_1_395	\N	strain	sample1_elastomer1	395	1.97	\N	0
mechanic_stress_sample1_elastomer1_tens_1_85	\N	mechanic_stress	sample1_elastomer1	85	356300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_114	\N	easy_hard_axis	Ni2MnGa_sample_1	114	0	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_1	\N	easy_hard_axis	Ni2MnGa_sample_1	1	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_20	\N	magnetic_field_strength	Ni2MnGa_sample_1	20	190000	\N	0
strain_sample1_elastomer1_tens_1_302	\N	strain	sample1_elastomer1	302	1.505	\N	0
strain_sample1_elastomer1_tens_1_199	\N	strain	sample1_elastomer1	199	0.99	\N	0
mechanic_stress_sample1_elastomer1_tens_1_339	\N	mechanic_stress	sample1_elastomer1	339	1088300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_7	\N	magnetic_field_strength	Ni2MnGa_sample_1	7	60000	\N	0
strain_sample1_elastomer1_tens_1_367	\N	strain	sample1_elastomer1	367	1.83	\N	0
strain_sample1_elastomer1_tens_1_228	\N	strain	sample1_elastomer1	228	1.135	\N	0
martensite_content_wire_actuator_NiTi#6_214	\N	martensite_content	wire_actuator_NiTi#6	214	0.5	\N	0
strain_sample1_elastomer1_tens_1_286	\N	strain	sample1_elastomer1	286	1.425	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_60	\N	magnetic_field_strength	Ni2MnGa_sample_1	60	590000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_43	\N	mechanic_stress	sample1_elastomer1	43	213300	\N	0
strain_sample1_elastomer1_tens_1_63	\N	strain	sample1_elastomer1	63	0.31	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_109	\N	magnetic_field_strength	Ni2MnGa_sample_1	109	80000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_317	\N	mechanic_stress	sample1_elastomer1	317	1008300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_14	\N	magnetic_field_strength	Ni2MnGa_sample_1	14	130000	\N	0
strain_sample1_elastomer1_tens_1_332	\N	strain	sample1_elastomer1	332	1.655	\N	0
mechanic_stress_sample1_elastomer1_tens_1_171	\N	mechanic_stress	sample1_elastomer1	171	586300	\N	0
strain_sample1_elastomer1_tens_1_323	\N	strain	sample1_elastomer1	323	1.61	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_52	\N	easy_hard_axis	Ni2MnGa_sample_1	52	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_12	\N	mechanic_stress	sample1_elastomer1	12	64900	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_129	\N	easy_hard_axis	Ni2MnGa_sample_1	129	0	\N	0
strain_sample1_elastomer1_tens_1_183	\N	strain	sample1_elastomer1	183	0.91	\N	0
hysteresis_wire_actuator_1_136	\N	hysteresis	wire_actuator_1	136	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_79	\N	mechanic_stress	sample1_elastomer1	79	339300	\N	0
dielectric_strength_M420	\N	dielectric_strength	M420	0	2000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_330	\N	mechanic_stress	sample1_elastomer1	330	1053300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_197	\N	mechanic_stress	sample1_elastomer1	197	654300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_3	\N	mechanic_stress	sample1_elastomer1	3	12150	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_36	\N	easy_hard_axis	Ni2MnGa_sample_1	36	1	\N	0
strain_sample1_elastomer1_tens_1_320	\N	strain	sample1_elastomer1	320	1.595	\N	0
strain_sample1_elastomer1_tens_1_84	\N	strain	sample1_elastomer1	84	0.415	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_133	\N	magnetic_field_strength	Ni2MnGa_sample_1	133	320000	\N	0
strain_sample1_elastomer1_tens_1_58	\N	strain	sample1_elastomer1	58	0.285	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_122	\N	easy_hard_axis	Ni2MnGa_sample_1	122	0	\N	0
strain_sample1_elastomer1_tens_1_49	\N	strain	sample1_elastomer1	49	0.24	\N	0
mechanic_stress_sample1_elastomer1_tens_1_78	\N	mechanic_stress	sample1_elastomer1	78	337300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_202	\N	mechanic_stress	sample1_elastomer1	202	666800	\N	0
mechanic_stress_sample1_elastomer1_tens_1_318	\N	mechanic_stress	sample1_elastomer1	318	1013300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_50	\N	mechanic_stress	sample1_elastomer1	50	241300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_219	\N	mechanic_stress	sample1_elastomer1	219	712300	\N	0
strain_sample1_elastomer1_tens_1_278	\N	strain	sample1_elastomer1	278	1.385	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_161	\N	magnetic_field_strength	Ni2MnGa_sample_1	161	600000	\N	0
strain_sample1_elastomer1_tens_1_352	\N	strain	sample1_elastomer1	352	1.755	\N	0
strain_sample1_elastomer1_tens_1_288	\N	strain	sample1_elastomer1	288	1.435	\N	0
strain_sample1_elastomer1_tens_1_173	\N	strain	sample1_elastomer1	173	0.86	\N	0
mechanic_stress_sample1_elastomer1_tens_1_59	\N	mechanic_stress	sample1_elastomer1	59	273300	\N	0
strain_sample1_elastomer1_tens_1_312	\N	strain	sample1_elastomer1	312	1.555	\N	0
strain_sample1_elastomer1_tens_1_353	\N	strain	sample1_elastomer1	353	1.76	\N	0
strain_sample1_elastomer1_tens_1_221	\N	strain	sample1_elastomer1	221	1.1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_247	\N	mechanic_stress	sample1_elastomer1	247	792300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_174	\N	mechanic_stress	sample1_elastomer1	174	594300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_4	\N	easy_hard_axis	Ni2MnGa_sample_1	4	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_315	\N	mechanic_stress	sample1_elastomer1	315	1008300	\N	0
\N	\N	vacuum_permittivity	sample1_elastomer1	0	8.854e-12	\N	0
martensite_content_wire_actuator_NiTi#6_316	\N	martensite_content	wire_actuator_NiTi#6	316	1	\N	0
martensite_content_wire_actuator_NiTi#6_221	\N	martensite_content	wire_actuator_NiTi#6	221	0.5	\N	0
strain_sample1_elastomer1_tens_1_55	\N	strain	sample1_elastomer1	55	0.27	\N	0
strain_sample1_elastomer1_tens_1_351	\N	strain	sample1_elastomer1	351	1.75	\N	0
strain_sample1_elastomer1_tens_1_209	\N	strain	sample1_elastomer1	209	1.04	\N	0
mechanic_stress_sample1_elastomer1_tens_1_289	\N	mechanic_stress	sample1_elastomer1	289	920300	\N	0
strain_sample1_elastomer1_tens_1_137	\N	strain	sample1_elastomer1	137	0.68	\N	0
austenite_start_temperature_wire_actuator_1_0	\N	austenite_start_temperature	wire_actuator_1	0	69.4	\N	0
strain_sample1_elastomer1_tens_1_126	\N	strain	sample1_elastomer1	126	0.625	\N	0
mechanic_stress_sample1_elastomer1_tens_1_348	\N	mechanic_stress	sample1_elastomer1	348	1118300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_11	\N	magnetic_field_strength	Ni2MnGa_sample_1	11	100000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_108	\N	magnetic_field_strength	Ni2MnGa_sample_1	108	70000	\N	0
strain_sample1_elastomer1_tens_1_339	\N	strain	sample1_elastomer1	339	1.69	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_2	\N	magnetic_field_strength	Ni2MnGa_sample_1	2	10000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_139	\N	mechanic_stress	sample1_elastomer1	139	503300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_60	\N	mechanic_stress	sample1_elastomer1	60	277300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_311	\N	mechanic_stress	sample1_elastomer1	311	988300	\N	0
temperature_wire_actuator_1_131	\N	temperature	wire_actuator_1	131	45	\N	0
mechanic_stress_sample1_elastomer1_tens_1_192	\N	mechanic_stress	sample1_elastomer1	192	641300	\N	0
strain_sample1_elastomer1_tens_1_45	\N	strain	sample1_elastomer1	45	0.22	\N	0
mechanic_stress_sample1_elastomer1_tens_1_298	\N	mechanic_stress	sample1_elastomer1	298	948800	\N	0
strain_wire_actuator_NiTi#6_311	\N	strain	wire_actuator_NiTi#6	311	0.005	\N	0
strain_sample1_elastomer1_tens_1_226	\N	strain	sample1_elastomer1	226	1.125	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_44	\N	magnetic_field_strength	Ni2MnGa_sample_1	44	430000	\N	0
martensite_content_wire_actuator_NiTi#6_310	\N	martensite_content	wire_actuator_NiTi#6	310	1	\N	0
strain_wire_actuator_NiTi#6_322	\N	strain	wire_actuator_NiTi#6	322	0.06	\N	0
strain_wire_actuator_NiTi#6_122	\N	strain	wire_actuator_NiTi#6	122	0.06	\N	0
strain_sample1_elastomer1_tens_1_28	\N	strain	sample1_elastomer1	28	0.135	\N	0
strain_wire_actuator_NiTi#6_212	\N	strain	wire_actuator_NiTi#6	212	0.01	\N	0
strain_wire_actuator_NiTi#6_318	\N	strain	wire_actuator_NiTi#6	318	0.04	\N	0
mechanic_stress_sample1_elastomer1_tens_1_319	\N	mechanic_stress	sample1_elastomer1	319	1018300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_32	\N	magnetic_field_strength	Ni2MnGa_sample_1	32	310000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_250	\N	mechanic_stress	sample1_elastomer1	250	800800	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_130	\N	magnetic_field_strength	Ni2MnGa_sample_1	130	290000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_121	\N	magnetic_field_strength	Ni2MnGa_sample_1	121	200000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_342	\N	mechanic_stress	sample1_elastomer1	342	1093300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_139	\N	magnetic_field_strength	Ni2MnGa_sample_1	139	380000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_138	\N	easy_hard_axis	Ni2MnGa_sample_1	138	0	\N	0
strain_sample1_elastomer1_tens_1_125	\N	strain	sample1_elastomer1	125	0.62	\N	0
strain_sample1_elastomer1_tens_1_115	\N	strain	sample1_elastomer1	115	0.57	\N	0
strain_sample1_elastomer1_tens_1_285	\N	strain	sample1_elastomer1	285	1.42	\N	0
strain_sample1_elastomer1_tens_1_208	\N	strain	sample1_elastomer1	208	1.035	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_29	\N	magnetic_field_strength	Ni2MnGa_sample_1	29	280000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_303	\N	mechanic_stress	sample1_elastomer1	303	965300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_121	\N	mechanic_stress	sample1_elastomer1	121	457300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_75	\N	mechanic_stress	sample1_elastomer1	75	327300	\N	0
strain_sample1_elastomer1_tens_1_54	\N	strain	sample1_elastomer1	54	0.265	\N	0
strain_sample1_elastomer1_tens_1_217	\N	strain	sample1_elastomer1	217	1.08	\N	0
strain_sample1_elastomer1_tens_1_380	\N	strain	sample1_elastomer1	380	1.895	\N	0
temperature_wire_actuator_1_135	\N	temperature	wire_actuator_1	135	52	\N	0
mechanic_stress_sample1_elastomer1_tens_1_384	\N	mechanic_stress	sample1_elastomer1	384	1243300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_56	\N	easy_hard_axis	Ni2MnGa_sample_1	56	1	\N	0
strain_sample1_elastomer1_tens_1_37	\N	strain	sample1_elastomer1	37	0.18	\N	0
mechanic_stress_sample1_elastomer1_tens_1_325	\N	mechanic_stress	sample1_elastomer1	325	1038300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_28	\N	easy_hard_axis	Ni2MnGa_sample_1	28	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_184	\N	mechanic_stress	sample1_elastomer1	184	621300	\N	0
strain_sample1_elastomer1_tens_1_92	\N	strain	sample1_elastomer1	92	0.455	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_118	\N	magnetic_field_strength	Ni2MnGa_sample_1	118	170000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_37	\N	mechanic_stress	sample1_elastomer1	37	187300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_18	\N	easy_hard_axis	Ni2MnGa_sample_1	18	1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_146	\N	magnetic_field_strength	Ni2MnGa_sample_1	146	450000	\N	0
M420_1_6	\N	B00	M420	0	-1.6e-10	\N	0
strain_sample1_elastomer1_tens_1_200	\N	strain	sample1_elastomer1	200	0.995	\N	0
mechanic_stress_sample1_elastomer1_tens_1_267	\N	mechanic_stress	sample1_elastomer1	267	851300	\N	0
strain_sample1_elastomer1_tens_1_256	\N	strain	sample1_elastomer1	256	1.275	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_51	\N	easy_hard_axis	Ni2MnGa_sample_1	51	1	\N	0
hysteresis_wire_actuator_1_134	\N	hysteresis	wire_actuator_1	134	1	\N	0
strain_sample1_elastomer1_tens_1_184	\N	strain	sample1_elastomer1	184	0.915	\N	0
strain_sample1_elastomer1_tens_1_248	\N	strain	sample1_elastomer1	248	1.235	\N	0
hysteresis_wire_actuator_1_118	\N	hysteresis	wire_actuator_1	118	2	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_103	\N	easy_hard_axis	Ni2MnGa_sample_1	103	0	\N	0
M420_1_3	\N	A00	M420	0	-6.5e-12	\N	0
martensite_content_wire_actuator_NiTi#6_110	\N	martensite_content	wire_actuator_NiTi#6	110	0	\N	0
temperature_wire_actuator_1_114	\N	temperature	wire_actuator_1	114	51	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_139	\N	easy_hard_axis	Ni2MnGa_sample_1	139	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_238	\N	mechanic_stress	sample1_elastomer1	238	765800	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_9	\N	easy_hard_axis	Ni2MnGa_sample_1	9	1	\N	0
martensite_content_wire_actuator_NiTi#6_210	\N	martensite_content	wire_actuator_NiTi#6	210	0.5	\N	0
strain_sample1_elastomer1_tens_1_1	\N	strain	sample1_elastomer1	1	0	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_18	\N	magnetic_field_strength	Ni2MnGa_sample_1	18	170000	\N	0
strain_sample1_elastomer1_tens_1_399	\N	strain	sample1_elastomer1	399	1.99	\N	0
temperature_wire_actuator_1_137	\N	temperature	wire_actuator_1	137	60	\N	0
strain_sample1_elastomer1_tens_1_48	\N	strain	sample1_elastomer1	48	0.235	\N	0
strain_sample1_elastomer1_tens_1_333	\N	strain	sample1_elastomer1	333	1.66	\N	0
strain_sample1_elastomer1_tens_1_232	\N	strain	sample1_elastomer1	232	1.155	\N	0
mechanic_stress_sample1_elastomer1_tens_1_128	\N	mechanic_stress	sample1_elastomer1	128	475300	\N	0
plateau_finish_strain_wire_actuator_NiTi#6_0	\N	plateau_finish_strain	wire_actuator_NiTi#6	0	0.0581	\N	0
strain_wire_actuator_NiTi#6_313	\N	strain	wire_actuator_NiTi#6	313	0.015	\N	0
mechanic_stress_sample1_elastomer1_tens_1_187	\N	mechanic_stress	sample1_elastomer1	187	628300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_349	\N	mechanic_stress	sample1_elastomer1	349	1118300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_164	\N	mechanic_stress	sample1_elastomer1	164	568300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_141	\N	easy_hard_axis	Ni2MnGa_sample_1	141	0	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_17	\N	easy_hard_axis	Ni2MnGa_sample_1	17	1	\N	0
martensite_content_wire_actuator_NiTi#6_218	\N	martensite_content	wire_actuator_NiTi#6	218	0.5	\N	0
mechanic_stress_sample1_elastomer1_tens_1_162	\N	mechanic_stress	sample1_elastomer1	162	562300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_141	\N	mechanic_stress	sample1_elastomer1	141	509300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_158	\N	magnetic_field_strength	Ni2MnGa_sample_1	158	570000	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_28	\N	magnetic_field_strength	Ni2MnGa_sample_1	28	270000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_292	\N	mechanic_stress	sample1_elastomer1	292	930300	\N	0
strain_sample1_elastomer1_tens_1_393	\N	strain	sample1_elastomer1	393	1.96	\N	0
mechanic_stress_sample1_elastomer1_tens_1_293	\N	mechanic_stress	sample1_elastomer1	293	933300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_29	\N	mechanic_stress	sample1_elastomer1	29	151300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_46	\N	magnetic_field_strength	Ni2MnGa_sample_1	46	450000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_32	\N	easy_hard_axis	Ni2MnGa_sample_1	32	1	\N	0
strain_sample1_elastomer1_tens_1_269	\N	strain	sample1_elastomer1	269	1.34	\N	0
strain_sample1_elastomer1_tens_1_43	\N	strain	sample1_elastomer1	43	0.21	\N	0
max_blocking_stress_NiTi	\N	maximal_blocking_stress	NiTi	0	500000000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_136	\N	mechanic_stress	sample1_elastomer1	136	497300	\N	0
strain_sample1_elastomer1_tens_1_180	\N	strain	sample1_elastomer1	180	0.895	\N	0
mechanic_stress_sample1_elastomer1_tens_1_266	\N	mechanic_stress	sample1_elastomer1	266	848300	\N	0
strain_wire_actuator_NiTi#6_315	\N	strain	wire_actuator_NiTi#6	315	0.025	\N	0
strain_sample1_elastomer1_tens_1_29	\N	strain	sample1_elastomer1	29	0.14	\N	0
mechanic_stress_sample1_elastomer1_tens_1_308	\N	mechanic_stress	sample1_elastomer1	308	983300	\N	0
strain_sample1_elastomer1_tens_1_370	\N	strain	sample1_elastomer1	370	1.845	\N	0
mechanic_stress_sample1_elastomer1_tens_1_395	\N	mechanic_stress	sample1_elastomer1	395	1278300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_324	\N	mechanic_stress	sample1_elastomer1	324	1033300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_185	\N	mechanic_stress	sample1_elastomer1	185	622300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_153	\N	magnetic_field_strength	Ni2MnGa_sample_1	153	520000	\N	0
hysteresis_wire_actuator_1_132	\N	hysteresis	wire_actuator_1	132	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_158	\N	mechanic_stress	sample1_elastomer1	158	553300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_154	\N	mechanic_stress	sample1_elastomer1	154	542300	\N	0
E-modulus_de-twinned_martensite_wire_actuator_NiTi#6_0	\N	E-modulus_de-twinned_martensite	wire_actuator_NiTi#6	0	24830000000	\N	0
strain_wire_actuator_NiTi#6_118	\N	strain	wire_actuator_NiTi#6	118	0.04	\N	0
strain_sample1_elastomer1_tens_1_78	\N	strain	sample1_elastomer1	78	0.385	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_124	\N	magnetic_field_strength	Ni2MnGa_sample_1	124	230000	\N	0
strain_sample1_elastomer1_tens_1_298	\N	strain	sample1_elastomer1	298	1.485	\N	0
strain_sample1_elastomer1_tens_1_306	\N	strain	sample1_elastomer1	306	1.525	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_14	\N	easy_hard_axis	Ni2MnGa_sample_1	14	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_396	\N	mechanic_stress	sample1_elastomer1	396	1283300	\N	0
strain_sample1_elastomer1_tens_1_250	\N	strain	sample1_elastomer1	250	1.245	\N	0
strain_sample1_elastomer1_tens_1_238	\N	strain	sample1_elastomer1	238	1.185	\N	0
mechanic_stress_sample1_elastomer1_tens_1_288	\N	mechanic_stress	sample1_elastomer1	288	916800	\N	0
mechanic_stress_sample1_elastomer1_tens_1_240	\N	mechanic_stress	sample1_elastomer1	240	771300	\N	0
strain_sample1_elastomer1_tens_1_74	\N	strain	sample1_elastomer1	74	0.365	\N	0
martensite_content_wire_actuator_NiTi#6_317	\N	martensite_content	wire_actuator_NiTi#6	317	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_36	\N	mechanic_stress	sample1_elastomer1	36	183300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_37	\N	easy_hard_axis	Ni2MnGa_sample_1	37	1	\N	0
martensite_content_wire_actuator_NiTi#6_113	\N	martensite_content	wire_actuator_NiTi#6	113	0	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_130	\N	easy_hard_axis	Ni2MnGa_sample_1	130	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_296	\N	mechanic_stress	sample1_elastomer1	296	942300	\N	0
strain_sample1_elastomer1_tens_1_46	\N	strain	sample1_elastomer1	46	0.225	\N	0
mechanic_stress_sample1_elastomer1_tens_1_214	\N	mechanic_stress	sample1_elastomer1	214	698800	\N	0
strain_sample1_elastomer1_tens_1_309	\N	strain	sample1_elastomer1	309	1.54	\N	0
strain_sample1_elastomer1_tens_1_390	\N	strain	sample1_elastomer1	390	1.945	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_61	\N	magnetic_field_strength	Ni2MnGa_sample_1	61	600000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_133	\N	mechanic_stress	sample1_elastomer1	133	488300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_109	\N	easy_hard_axis	Ni2MnGa_sample_1	109	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_216	\N	mechanic_stress	sample1_elastomer1	216	704300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_66	\N	mechanic_stress	sample1_elastomer1	66	297300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_30	\N	mechanic_stress	sample1_elastomer1	30	156300	\N	0
temperature_wire_actuator_1_111	\N	temperature	wire_actuator_1	111	45	\N	0
martensite_content_wire_actuator_NiTi#6_315	\N	martensite_content	wire_actuator_NiTi#6	315	1	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_112	\N	easy_hard_axis	Ni2MnGa_sample_1	112	0	\N	0
strain_sample1_elastomer1_tens_1_117	\N	strain	sample1_elastomer1	117	0.58	\N	0
strain_sample1_elastomer1_tens_1_215	\N	strain	sample1_elastomer1	215	1.07	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_21	\N	magnetic_field_strength	Ni2MnGa_sample_1	21	200000	\N	0
strain_sample1_elastomer1_tens_1_97	\N	strain	sample1_elastomer1	97	0.48	\N	0
strain_wire_actuator_NiTi#6_116	\N	strain	wire_actuator_NiTi#6	116	0.03	\N	0
mechanic_stress_sample1_elastomer1_tens_1_120	\N	mechanic_stress	sample1_elastomer1	120	454300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_314	\N	mechanic_stress	sample1_elastomer1	314	1003300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_80	\N	mechanic_stress	sample1_elastomer1	80	342300	\N	0
strain_sample1_elastomer1_tens_1_155	\N	strain	sample1_elastomer1	155	0.77	\N	0
mechanic_stress_sample1_elastomer1_tens_1_182	\N	mechanic_stress	sample1_elastomer1	182	613300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_389	\N	mechanic_stress	sample1_elastomer1	389	1258300	\N	0
strain_sample1_elastomer1_tens_1_230	\N	strain	sample1_elastomer1	230	1.145	\N	0
mechanic_stress_sample1_elastomer1_tens_1_284	\N	mechanic_stress	sample1_elastomer1	284	904300	\N	0
strain_sample1_elastomer1_tens_1_143	\N	strain	sample1_elastomer1	143	0.71	\N	0
mechanic_stress_sample1_elastomer1_tens_1_124	\N	mechanic_stress	sample1_elastomer1	124	464300	\N	0
strain_sample1_elastomer1_tens_1_374	\N	strain	sample1_elastomer1	374	1.865	\N	0
mechanic_stress_sample1_elastomer1_tens_1_101	\N	mechanic_stress	sample1_elastomer1	101	402300	\N	0
strain_sample1_elastomer1_tens_1_216	\N	strain	sample1_elastomer1	216	1.075	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_156	\N	easy_hard_axis	Ni2MnGa_sample_1	156	0	\N	0
strain_sample1_elastomer1_tens_1_335	\N	strain	sample1_elastomer1	335	1.67	\N	0
strain_sample1_elastomer1_tens_1_240	\N	strain	sample1_elastomer1	240	1.195	\N	0
mechanic_stress_sample1_elastomer1_tens_1_245	\N	mechanic_stress	sample1_elastomer1	245	786300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_119	\N	magnetic_field_strength	Ni2MnGa_sample_1	119	180000	\N	0
strain_sample1_elastomer1_tens_1_161	\N	strain	sample1_elastomer1	161	0.8	\N	0
strain_sample1_elastomer1_tens_1_142	\N	strain	sample1_elastomer1	142	0.705	\N	0
mechanic_stress_sample1_elastomer1_tens_1_127	\N	mechanic_stress	sample1_elastomer1	127	474300	\N	0
hysteresis_wire_actuator_1_138	\N	hysteresis	wire_actuator_1	138	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_167	\N	mechanic_stress	sample1_elastomer1	167	576300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_34	\N	magnetic_field_strength	Ni2MnGa_sample_1	34	330000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_151	\N	easy_hard_axis	Ni2MnGa_sample_1	151	0	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_60	\N	easy_hard_axis	Ni2MnGa_sample_1	60	1	\N	0
strain_sample1_elastomer1_tens_1_198	\N	strain	sample1_elastomer1	198	0.985	\N	0
strain_sample1_elastomer1_tens_1_210	\N	strain	sample1_elastomer1	210	1.045	\N	0
strain_sample1_elastomer1_tens_1_96	\N	strain	sample1_elastomer1	96	0.475	\N	0
strain_sample1_elastomer1_tens_1_283	\N	strain	sample1_elastomer1	283	1.41	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_35	\N	magnetic_field_strength	Ni2MnGa_sample_1	35	340000	\N	0
strain_sample1_elastomer1_tens_1_147	\N	strain	sample1_elastomer1	147	0.73	\N	0
strain_sample1_elastomer1_tens_1_360	\N	strain	sample1_elastomer1	360	1.795	\N	0
strain_sample1_elastomer1_tens_1_64	\N	strain	sample1_elastomer1	64	0.315	\N	0
strain_sample1_elastomer1_tens_1_222	\N	strain	sample1_elastomer1	222	1.105	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_143	\N	magnetic_field_strength	Ni2MnGa_sample_1	143	420000	\N	0
\N	\N	vacuum_permittivity	M420	0	8.854e-12	\N	0
\N	\N	vacuum_permeability	NiTi	0	1.2566e-06	\N	0
strain_sample1_elastomer1_tens_1_244	\N	strain	sample1_elastomer1	244	1.215	\N	0
mechanic_stress_sample1_elastomer1_tens_1_116	\N	mechanic_stress	sample1_elastomer1	116	443300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_243	\N	mechanic_stress	sample1_elastomer1	243	778300	\N	0
martensite_content_wire_actuator_NiTi#6_319	\N	martensite_content	wire_actuator_NiTi#6	319	1	\N	0
strain_sample1_elastomer1_tens_1_342	\N	strain	sample1_elastomer1	342	1.705	\N	0
martensite_content_wire_actuator_NiTi#6_313	\N	martensite_content	wire_actuator_NiTi#6	313	1	\N	0
strain_sample1_elastomer1_tens_1_138	\N	strain	sample1_elastomer1	138	0.685	\N	0
mechanic_stress_sample1_elastomer1_tens_1_220	\N	mechanic_stress	sample1_elastomer1	220	715300	\N	0
hysteresis_wire_actuator_1_133	\N	hysteresis	wire_actuator_1	133	1	\N	0
strain_sample1_elastomer1_tens_1_384	\N	strain	sample1_elastomer1	384	1.915	\N	0
strain_sample1_elastomer1_tens_1_33	\N	strain	sample1_elastomer1	33	0.16	\N	0
hysteresis_wire_actuator_1_131	\N	hysteresis	wire_actuator_1	131	1	\N	0
strain_sample1_elastomer1_tens_1_178	\N	strain	sample1_elastomer1	178	0.885	\N	0
strain_sample1_elastomer1_tens_1_13	\N	strain	sample1_elastomer1	13	0.06	\N	0
mechanic_stress_sample1_elastomer1_tens_1_353	\N	mechanic_stress	sample1_elastomer1	353	1128300	\N	0
strain_sample1_elastomer1_tens_1_169	\N	strain	sample1_elastomer1	169	0.84	\N	0
hysteresis_wire_actuator_1_116	\N	hysteresis	wire_actuator_1	116	2	\N	0
strain_sample1_elastomer1_tens_1_378	\N	strain	sample1_elastomer1	378	1.885	\N	0
mechanic_stress_sample1_elastomer1_tens_1_100	\N	mechanic_stress	sample1_elastomer1	100	400300	\N	0
strain_sample1_elastomer1_tens_1_150	\N	strain	sample1_elastomer1	150	0.745	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_161	\N	easy_hard_axis	Ni2MnGa_sample_1	161	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_256	\N	mechanic_stress	sample1_elastomer1	256	818300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_55	\N	magnetic_field_strength	Ni2MnGa_sample_1	55	540000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_64	\N	mechanic_stress	sample1_elastomer1	64	290300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_155	\N	easy_hard_axis	Ni2MnGa_sample_1	155	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_394	\N	mechanic_stress	sample1_elastomer1	394	1273300	\N	0
strain_sample1_elastomer1_tens_1_187	\N	strain	sample1_elastomer1	187	0.93	\N	0
strain_sample1_elastomer1_tens_1_22	\N	strain	sample1_elastomer1	22	0.105	\N	0
strain_sample1_elastomer1_tens_1_227	\N	strain	sample1_elastomer1	227	1.13	\N	0
strain_sample1_elastomer1_tens_1_130	\N	strain	sample1_elastomer1	130	0.645	\N	0
mechanic_stress_sample1_elastomer1_tens_1_209	\N	mechanic_stress	sample1_elastomer1	209	685300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_391	\N	mechanic_stress	sample1_elastomer1	391	1268300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_225	\N	mechanic_stress	sample1_elastomer1	225	729300	\N	0
martensite_start_temperature_wire_actuator_1_0	\N	martensite_start_temperature	wire_actuator_1	0	52.4	\N	0
strain_wire_actuator_NiTi#6_317	\N	strain	wire_actuator_NiTi#6	317	0.035	\N	0
strain_wire_actuator_NiTi#6_114	\N	strain	wire_actuator_NiTi#6	114	0.02	\N	0
strain_sample1_elastomer1_tens_1_311	\N	strain	sample1_elastomer1	311	1.55	\N	0
strain_sample1_elastomer1_tens_1_249	\N	strain	sample1_elastomer1	249	1.24	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_111	\N	magnetic_field_strength	Ni2MnGa_sample_1	111	100000	\N	0
strain_sample1_elastomer1_tens_1_24	\N	strain	sample1_elastomer1	24	0.115	\N	0
strain_sample1_elastomer1_tens_1_355	\N	strain	sample1_elastomer1	355	1.77	\N	0
mechanic_stress_sample1_elastomer1_tens_1_188	\N	mechanic_stress	sample1_elastomer1	188	630300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_376	\N	mechanic_stress	sample1_elastomer1	376	1213300	\N	0
hysteresis_wire_actuator_1_123	\N	hysteresis	wire_actuator_1	123	2	\N	0
mechanic_stress_sample1_elastomer1_tens_1_114	\N	mechanic_stress	sample1_elastomer1	114	438300	\N	0
strain_wire_actuator_NiTi#6_112	\N	strain	wire_actuator_NiTi#6	112	0.01	\N	0
mechanic_stress_sample1_elastomer1_tens_1_199	\N	mechanic_stress	sample1_elastomer1	199	659300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_10	\N	easy_hard_axis	Ni2MnGa_sample_1	10	1	\N	0
mechanic_stress_sample1_elastomer1_tens_1_91	\N	mechanic_stress	sample1_elastomer1	91	375300	\N	0
strain_sample1_elastomer1_tens_1_337	\N	strain	sample1_elastomer1	337	1.68	\N	0
strain_sample1_elastomer1_tens_1_372	\N	strain	sample1_elastomer1	372	1.855	\N	0
mechanic_stress_sample1_elastomer1_tens_1_251	\N	mechanic_stress	sample1_elastomer1	251	803300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_221	\N	mechanic_stress	sample1_elastomer1	221	718300	\N	0
martensite_content_wire_actuator_NiTi#6_120	\N	martensite_content	wire_actuator_NiTi#6	120	0	\N	0
strain_sample1_elastomer1_tens_1_233	\N	strain	sample1_elastomer1	233	1.16	\N	0
strain_sample1_elastomer1_tens_1_219	\N	strain	sample1_elastomer1	219	1.09	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_145	\N	easy_hard_axis	Ni2MnGa_sample_1	145	0	\N	0
strain_sample1_elastomer1_tens_1_182	\N	strain	sample1_elastomer1	182	0.905	\N	0
mechanic_stress_sample1_elastomer1_tens_1_304	\N	mechanic_stress	sample1_elastomer1	304	968800	\N	0
mechanic_stress_sample1_elastomer1_tens_1_273	\N	mechanic_stress	sample1_elastomer1	273	871300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_155	\N	magnetic_field_strength	Ni2MnGa_sample_1	155	540000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_38	\N	easy_hard_axis	Ni2MnGa_sample_1	38	1	\N	0
hysteresis_wire_actuator_1_143	\N	hysteresis	wire_actuator_1	143	1	\N	0
strain_sample1_elastomer1_tens_1_124	\N	strain	sample1_elastomer1	124	0.615	\N	0
mechanic_stress_sample1_elastomer1_tens_1_370	\N	mechanic_stress	sample1_elastomer1	370	1193300	\N	0
strain_sample1_elastomer1_tens_1_99	\N	strain	sample1_elastomer1	99	0.49	\N	0
strain_sample1_elastomer1_tens_1_112	\N	strain	sample1_elastomer1	112	0.555	\N	0
strain_wire_actuator_NiTi#6_215	\N	strain	wire_actuator_NiTi#6	215	0.025	\N	0
E-modulus_austenite_wire_actuator_NiTi#6_0	\N	E-modulus_austenite	wire_actuator_NiTi#6	0	53510000000	\N	0
strain_sample1_elastomer1_tens_1_21	\N	strain	sample1_elastomer1	21	0.1	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_48	\N	magnetic_field_strength	Ni2MnGa_sample_1	48	470000	\N	0
strain_sample1_elastomer1_tens_1_108	\N	strain	sample1_elastomer1	108	0.535	\N	0
mechanic_stress_sample1_elastomer1_tens_1_65	\N	mechanic_stress	sample1_elastomer1	65	294300	\N	0
strain_wire_actuator_NiTi#6_218	\N	strain	wire_actuator_NiTi#6	218	0.04	\N	0
mechanic_stress_sample1_elastomer1_tens_1_18	\N	mechanic_stress	sample1_elastomer1	18	97967	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_7	\N	easy_hard_axis	Ni2MnGa_sample_1	7	1	\N	0
\N	\N	vacuum_permeability	M420	0	1.2566e-06	\N	0
strain_sample1_elastomer1_tens_1_181	\N	strain	sample1_elastomer1	181	0.9	\N	0
strain_sample1_elastomer1_tens_1_290	\N	strain	sample1_elastomer1	290	1.445	\N	0
mechanic_stress_sample1_elastomer1_tens_1_11	\N	mechanic_stress	sample1_elastomer1	11	59400	\N	0
strain_sample1_elastomer1_tens_1_330	\N	strain	sample1_elastomer1	330	1.645	\N	0
mechanic_stress_sample1_elastomer1_tens_1_399	\N	mechanic_stress	sample1_elastomer1	399	1298300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_222	\N	mechanic_stress	sample1_elastomer1	222	720800	\N	0
strain_sample1_elastomer1_tens_1_266	\N	strain	sample1_elastomer1	266	1.325	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_10	\N	magnetic_field_strength	Ni2MnGa_sample_1	10	90000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_88	\N	mechanic_stress	sample1_elastomer1	88	366300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_183	\N	mechanic_stress	sample1_elastomer1	183	618300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_104	\N	easy_hard_axis	Ni2MnGa_sample_1	104	0	\N	0
strain_sample1_elastomer1_tens_1_196	\N	strain	sample1_elastomer1	196	0.975	\N	0
mechanic_stress_sample1_elastomer1_tens_1_110	\N	mechanic_stress	sample1_elastomer1	110	428300	\N	0
strain_sample1_elastomer1_tens_1_359	\N	strain	sample1_elastomer1	359	1.79	\N	0
mechanic_stress_sample1_elastomer1_tens_1_383	\N	mechanic_stress	sample1_elastomer1	383	1238300	\N	0
strain_sample1_elastomer1_tens_1_236	\N	strain	sample1_elastomer1	236	1.175	\N	0
strain_sample1_elastomer1_tens_1_7	\N	strain	sample1_elastomer1	7	0.03	\N	0
mechanic_stress_sample1_elastomer1_tens_1_320	\N	mechanic_stress	sample1_elastomer1	320	1023300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_357	\N	mechanic_stress	sample1_elastomer1	357	1148300	\N	0
strain_sample1_elastomer1_tens_1_87	\N	strain	sample1_elastomer1	87	0.43	\N	0
strain_sample1_elastomer1_tens_1_12	\N	strain	sample1_elastomer1	12	0.055	\N	0
strain_sample1_elastomer1_tens_1_241	\N	strain	sample1_elastomer1	241	1.2	\N	0
strain_sample1_elastomer1_tens_1_231	\N	strain	sample1_elastomer1	231	1.15	\N	0
mechanic_stress_sample1_elastomer1_tens_1_20	\N	mechanic_stress	sample1_elastomer1	20	108633	\N	0
\N	\N	vacuum_permittivity	elastomer_1	0	8.854e-12	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_62	\N	magnetic_field_strength	Ni2MnGa_sample_1	62	610000	\N	0
strain_sample1_elastomer1_tens_1_44	\N	strain	sample1_elastomer1	44	0.215	\N	0
strain_sample1_elastomer1_tens_1_307	\N	strain	sample1_elastomer1	307	1.53	\N	0
strain_sample1_elastomer1_tens_1_70	\N	strain	sample1_elastomer1	70	0.345	\N	0
strain_sample1_elastomer1_tens_1_41	\N	strain	sample1_elastomer1	41	0.2	\N	0
martensite_content_wire_actuator_NiTi#6_217	\N	martensite_content	wire_actuator_NiTi#6	217	0.5	\N	0
martensite_content_wire_actuator_NiTi#6_216	\N	martensite_content	wire_actuator_NiTi#6	216	0.5	\N	0
strain_sample1_elastomer1_tens_1_354	\N	strain	sample1_elastomer1	354	1.765	\N	0
mechanic_stress_sample1_elastomer1_tens_1_87	\N	mechanic_stress	sample1_elastomer1	87	362300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_69	\N	mechanic_stress	sample1_elastomer1	69	307300	\N	0
mechanic_stress_sample1_elastomer1_tens_1_329	\N	mechanic_stress	sample1_elastomer1	329	1048300	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_2	\N	easy_hard_axis	Ni2MnGa_sample_1	2	1	\N	0
strain_sample1_elastomer1_tens_1_386	\N	strain	sample1_elastomer1	386	1.925	\N	0
mechanic_stress_sample1_elastomer1_tens_1_354	\N	mechanic_stress	sample1_elastomer1	354	1133300	\N	0
strain_sample1_elastomer1_tens_1_82	\N	strain	sample1_elastomer1	82	0.405	\N	0
strain_sample1_elastomer1_tens_1_365	\N	strain	sample1_elastomer1	365	1.82	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_1_45	\N	magnetic_field_strength	Ni2MnGa_sample_1	45	440000	\N	0
mechanic_stress_sample1_elastomer1_tens_1_279	\N	mechanic_stress	sample1_elastomer1	279	889300	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_103	\N	magnetic_field_strength	Ni2MnGa_sample_1	103	20000	\N	0
temperature_wire_actuator_1_123	\N	temperature	wire_actuator_1	123	75	\N	0
strain_sample1_elastomer1_tens_1_289	\N	strain	sample1_elastomer1	289	1.44	\N	0
mechanic_stress_sample1_elastomer1_tens_1_335	\N	mechanic_stress	sample1_elastomer1	335	1078300	\N	0
strain_sample1_elastomer1_tens_1_148	\N	strain	sample1_elastomer1	148	0.735	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_123	\N	easy_hard_axis	Ni2MnGa_sample_1	123	0	\N	0
strain_sample1_elastomer1_tens_1_129	\N	strain	sample1_elastomer1	129	0.64	\N	0
magnetic_field_strength_Ni2MnGa_sample1_meas_flux_over_field_2_132	\N	magnetic_field_strength	Ni2MnGa_sample_1	132	310000	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_117	\N	easy_hard_axis	Ni2MnGa_sample_1	117	0	\N	0
martensite_content_wire_actuator_NiTi#6_121	\N	martensite_content	wire_actuator_NiTi#6	121	0	\N	0
mechanic_stress_sample1_elastomer1_tens_1_237	\N	mechanic_stress	sample1_elastomer1	237	763300	\N	0
strain_sample1_elastomer1_tens_1_375	\N	strain	sample1_elastomer1	375	1.87	\N	0
mechanic_stress_sample1_elastomer1_tens_1_301	\N	mechanic_stress	sample1_elastomer1	301	959300	\N	0
M420_2_3	\N	B01	M420	0	0.037	\N	0
strain_sample1_elastomer1_tens_1_154	\N	strain	sample1_elastomer1	154	0.765	\N	0
strain_sample1_elastomer1_tens_1_118	\N	strain	sample1_elastomer1	118	0.585	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_2_154	\N	easy_hard_axis	Ni2MnGa_sample_1	154	0	\N	0
strain_sample1_elastomer1_tens_1_239	\N	strain	sample1_elastomer1	239	1.19	\N	0
mechanic_stress_sample1_elastomer1_tens_1_178	\N	mechanic_stress	sample1_elastomer1	178	604300	\N	0
strain_sample1_elastomer1_tens_1_177	\N	strain	sample1_elastomer1	177	0.88	\N	0
martensite_content_wire_actuator_NiTi#6_117	\N	martensite_content	wire_actuator_NiTi#6	117	0	\N	0
easy_hard_axis_Ni2MnGa_sample1_meas_flux_over_field_1_43	\N	easy_hard_axis	Ni2MnGa_sample_1	43	1	\N	0
strain_sample1_elastomer1_tens_1_203	\N	strain	sample1_elastomer1	203	1.01	\N	0
strain_sample1_elastomer1_tens_1_381	\N	strain	sample1_elastomer1	381	1.9	\N	0
strain_sample1_elastomer1_tens_1_83	\N	strain	sample1_elastomer1	83	0.41	\N	0
cp_mod_A01_matrix_model_4_M420_0_1	mod_matrix_model_4_A01_M420_0	A01	M420	0	1.35929071606054e-11	matrix_model_4	1
cp_mod_A01_matrix_model_4_M420_0_2	mod_matrix_model_4_A01_M420_0	A01	M420	0	-7.50709283939462e-12	matrix_model_4	1
cp_mod_A01_matrix_model_4_M420_0_3	mod_matrix_model_4_A01_M420_0	A01	M420	0	-2.49051276259318e-12	matrix_model_4	1
cp_mod_A01_matrix_model_4_M420_0_4	mod_matrix_model_4_A01_M420_0	A01	M420	0	9.80395019200361e-12	matrix_model_4	1
cp_mod_A01_matrix_model_4_M420_0_5	mod_matrix_model_4_A01_M420_0	A01	M420	0	2.55437514117913e-11	matrix_model_4	1
cp_mod_A10_matrix_model_1_M420_0_1	mod_matrix_model_1_A10_M420_0	A10	M420	0	120197869377.449	matrix_model_1	1
cp_mod_A10_matrix_model_1_M420_0_2	mod_matrix_model_1_A10_M420_0	A10	M420	0	72804504448.5392	matrix_model_1	1
cp_mod_A10_matrix_model_1_M420_0_3	mod_matrix_model_1_A10_M420_0	A10	M420	0	67086386624.0066	matrix_model_1	1
cp_mod_A10_matrix_model_1_M420_0_4	mod_matrix_model_1_A10_M420_0	A10	M420	0	100113530808.133	matrix_model_1	1
cp_mod_A10_matrix_model_1_M420_0_5	mod_matrix_model_1_A10_M420_0	A10	M420	0	22222222222.2222	matrix_model_1	1
cp_mod_B01_matrix_model_5_M420_0_3	mod_matrix_model_5_B01_M420_0	B01	M420	0	0.0370595211203976	matrix_model_5	1
cp_mod_B01_matrix_model_5_M420_0_1	mod_matrix_model_5_B01_M420_0	B01	M420	0	-0.0112943302462164	matrix_model_5	1
cp_mod_B01_matrix_model_5_M420_0_2	mod_matrix_model_5_B01_M420_0	B01	M420	0	0.0250592952337926	matrix_model_5	1
cp_mod_B10_matrix_model_2_M420_0_3	mod_matrix_model_2_B10_M420_0	B10	M420	0	11.6666666666667	matrix_model_2	1
cp_mod_B10_matrix_model_2_M420_0_1	mod_matrix_model_2_B10_M420_0	B10	M420	0	-7.06471256063577	matrix_model_2	1
cp_mod_B10_matrix_model_2_M420_0_2	mod_matrix_model_2_B10_M420_0	B10	M420	0	14.0726597172051	matrix_model_2	1
cp_mod_C01_matrix_model_6_M420_0_1	mod_matrix_model_6_C01_M420_0	C01	M420	0	70589564.0388525	matrix_model_6	1
cp_mod_C01_matrix_model_6_M420_0_2	mod_matrix_model_6_C01_M420_0	C01	M420	0	70589564.0388525	matrix_model_6	1
cp_mod_C10_matrix_model_3_M420_0_1	mod_matrix_model_3_C10_M420_0	C10	M420	0	8.0414e-09	matrix_model_3	1
cp_mod_C10_matrix_model_3_M420_0_2	mod_matrix_model_3_C10_M420_0	C10	M420	0	6.90989778098875e-09	matrix_model_3	1
cp_mod_linear_interpol_initial_stress_initial_stress_linear_interpolation_initial_stress_sample1_elastomer1_0	mod_linear_interpolation_initial_stress_initial_stress_sample1_elastomer1_0	initial_stress	sample1_elastomer1	0	59400	linear_interpolation_initial_stress	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_1_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_1_magnetization	magnetization	Ni2MnGa_sample_1	1	0	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_2_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_2_magnetization	magnetization	Ni2MnGa_sample_1	2	79792.73436256566	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_3_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_3_magnetization	magnetization	Ni2MnGa_sample_1	3	99308.04551965621	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_4_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_4_magnetization	magnetization	Ni2MnGa_sample_1	4	107944.8511857393	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_5_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_5_magnetization	magnetization	Ni2MnGa_sample_1	5	116581.64889384052	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_6_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_6_magnetization	magnetization	Ni2MnGa_sample_1	6	125218.45455992362	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_7_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_7_magnetization	magnetization	Ni2MnGa_sample_1	7	133855.2602260067	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_8_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_8_magnetization	magnetization	Ni2MnGa_sample_1	8	142492.07385007164	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_9_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_9_magnetization	magnetization	Ni2MnGa_sample_1	9	151128.87155817286	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_10_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_10_magnetization	magnetization	Ni2MnGa_sample_1	10	159765.68518223776	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_11_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_11_magnetization	magnetization	Ni2MnGa_sample_1	11	168402.47493235715	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_12_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_12_magnetization	magnetization	Ni2MnGa_sample_1	12	177039.2885564221	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_13_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_13_magnetization	magnetization	Ni2MnGa_sample_1	13	185676.1101384689	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_14_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_14_magnetization	magnetization	Ni2MnGa_sample_1	14	194312.90784657013	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_15_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_15_magnetization	magnetization	Ni2MnGa_sample_1	15	202949.70555467135	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_16_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_16_magnetization	magnetization	Ni2MnGa_sample_1	16	211586.50326277257	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_17_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_17_magnetization	magnetization	Ni2MnGa_sample_1	17	220223.32484481938	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_18_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_18_magnetization	magnetization	Ni2MnGa_sample_1	18	228860.13846888428	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_19_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_19_magnetization	magnetization	Ni2MnGa_sample_1	19	237496.8963870763	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_20_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_20_magnetization	magnetization	Ni2MnGa_sample_1	20	246133.72592710488	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_21_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_21_magnetization	magnetization	Ni2MnGa_sample_1	21	254770.5395511698	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_22_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_22_magnetization	magnetization	Ni2MnGa_sample_1	22	263407.3690911985	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_23_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_23_magnetization	magnetization	Ni2MnGa_sample_1	23	272044.1667992997	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_24_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_24_magnetization	magnetization	Ni2MnGa_sample_1	24	280680.9406334554	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_25_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_25_magnetization	magnetization	Ni2MnGa_sample_1	25	289317.7622155022	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_26_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_26_magnetization	magnetization	Ni2MnGa_sample_1	26	297954.55992360343	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_27_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_27_magnetization	magnetization	Ni2MnGa_sample_1	27	306591.38150565024	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_28_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_28_magnetization	magnetization	Ni2MnGa_sample_1	28	315228.20308769704	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_29_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_29_magnetization	magnetization	Ni2MnGa_sample_1	29	323864.9530479071	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_30_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_30_magnetization	magnetization	Ni2MnGa_sample_1	30	332501.7825879357	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_31_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_31_magnetization	magnetization	Ni2MnGa_sample_1	31	341138.5962120007	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_32_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_32_magnetization	magnetization	Ni2MnGa_sample_1	32	349775.4257520293	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_33_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_33_magnetization	magnetization	Ni2MnGa_sample_1	33	358412.18367022125	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_34_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_34_magnetization	magnetization	Ni2MnGa_sample_1	34	367049.02116823173	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_35_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_35_magnetization	magnetization	Ni2MnGa_sample_1	35	375685.81887633295	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_36_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_36_magnetization	magnetization	Ni2MnGa_sample_1	36	384322.6165844342	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_37_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_37_magnetization	magnetization	Ni2MnGa_sample_1	37	392959.4461244629	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_38_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_38_magnetization	magnetization	Ni2MnGa_sample_1	38	401596.2120006366	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_39_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_39_magnetization	magnetization	Ni2MnGa_sample_1	39	410233.00970873784	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_40_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_40_magnetization	magnetization	Ni2MnGa_sample_1	40	418869.80741683906	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_41_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_41_magnetization	magnetization	Ni2MnGa_sample_1	41	427506.6847047589	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_42_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_42_magnetization	magnetization	Ni2MnGa_sample_1	42	436143.48241286003	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_43_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_43_magnetization	magnetization	Ni2MnGa_sample_1	43	444780.28012096137	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_44_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_44_magnetization	magnetization	Ni2MnGa_sample_1	44	453417.0778290626	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_45_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_45_magnetization	magnetization	Ni2MnGa_sample_1	45	462053.7959573453	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_46_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_46_magnetization	magnetization	Ni2MnGa_sample_1	46	470690.67324526503	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_47_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_47_magnetization	magnetization	Ni2MnGa_sample_1	47	479327.47095336637	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_48_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_48_magnetization	magnetization	Ni2MnGa_sample_1	48	487964.2686614675	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_49_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_49_magnetization	magnetization	Ni2MnGa_sample_1	49	496601.0663695687	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_50_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_50_magnetization	magnetization	Ni2MnGa_sample_1	50	505237.94365748845	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_51_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_51_magnetization	magnetization	Ni2MnGa_sample_1	51	513874.7413655898	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_52_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_52_magnetization	magnetization	Ni2MnGa_sample_1	52	522511.539073691	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_53_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_53_magnetization	magnetization	Ni2MnGa_sample_1	53	531148.3367817921	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_54_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_54_magnetization	magnetization	Ni2MnGa_sample_1	54	539785.0549100749	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_55_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_55_magnetization	magnetization	Ni2MnGa_sample_1	55	548421.9321979946	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_56_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_56_magnetization	magnetization	Ni2MnGa_sample_1	56	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_57_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_57_magnetization	magnetization	Ni2MnGa_sample_1	57	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_58_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_58_magnetization	magnetization	Ni2MnGa_sample_1	58	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_59_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_59_magnetization	magnetization	Ni2MnGa_sample_1	59	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_60_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_60_magnetization	magnetization	Ni2MnGa_sample_1	60	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_61_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_61_magnetization	magnetization	Ni2MnGa_sample_1	61	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_62_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_62_magnetization	magnetization	Ni2MnGa_sample_1	62	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_101_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_101_magnetization	magnetization	Ni2MnGa_sample_1	101	0	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_102_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_102_magnetization	magnetization	Ni2MnGa_sample_1	102	336629.15804551967	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_103_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_103_magnetization	magnetization	Ni2MnGa_sample_1	103	476139.12939678505	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_104_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_104_magnetization	magnetization	Ni2MnGa_sample_1	104	477665.89208976604	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_105_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_105_magnetization	magnetization	Ni2MnGa_sample_1	105	479192.69457265636	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_106_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_106_magnetization	magnetization	Ni2MnGa_sample_1	106	480719.48113958305	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_107_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_107_magnetization	magnetization	Ni2MnGa_sample_1	107	482246.2597485278	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_108_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_108_magnetization	magnetization	Ni2MnGa_sample_1	108	483773.0622314181	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_109_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_109_magnetization	magnetization	Ni2MnGa_sample_1	109	485299.82492439914	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_110_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_110_magnetization	magnetization	Ni2MnGa_sample_1	110	486826.62740728946	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_111_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_111_magnetization	magnetization	Ni2MnGa_sample_1	111	488353.3901002706	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_112_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_112_magnetization	magnetization	Ni2MnGa_sample_1	112	489880.19258316094	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_113_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_113_magnetization	magnetization	Ni2MnGa_sample_1	113	491406.9711921057	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_114_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_114_magnetization	magnetization	Ni2MnGa_sample_1	114	492933.7577590323	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_115_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_115_magnetization	magnetization	Ni2MnGa_sample_1	115	494460.5681999045	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_116_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_116_magnetization	magnetization	Ni2MnGa_sample_1	116	495987.3229349038	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_117_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_117_magnetization	magnetization	Ni2MnGa_sample_1	117	497514.1254177941	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_118_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_118_magnetization	magnetization	Ni2MnGa_sample_1	118	499040.90402673883	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_119_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_119_magnetization	magnetization	Ni2MnGa_sample_1	119	500567.69059366547	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_120_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_120_magnetization	magnetization	Ni2MnGa_sample_1	120	502094.4612446284	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_121_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_121_magnetization	magnetization	Ni2MnGa_sample_1	121	503621.2796434824	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_122_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_122_magnetization	magnetization	Ni2MnGa_sample_1	122	505148.06621040904	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_123_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_123_magnetization	magnetization	Ni2MnGa_sample_1	123	506674.836861372	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_124_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_124_magnetization	magnetization	Ni2MnGa_sample_1	124	508201.6313862804	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_125_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_125_magnetization	magnetization	Ni2MnGa_sample_1	125	509728.3940792616	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_126_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_126_magnetization	magnetization	Ni2MnGa_sample_1	126	511255.2124781156	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_127_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_127_magnetization	magnetization	Ni2MnGa_sample_1	127	512781.95129715104	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_128_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_128_magnetization	magnetization	Ni2MnGa_sample_1	128	514308.76969600504	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_129_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_129_magnetization	magnetization	Ni2MnGa_sample_1	129	515835.58809485915	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_130_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_130_magnetization	magnetization	Ni2MnGa_sample_1	130	517362.3269138946	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_131_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_131_magnetization	magnetization	Ni2MnGa_sample_1	131	518889.0657329301	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_132_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_132_magnetization	magnetization	Ni2MnGa_sample_1	132	520415.96371160285	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_133_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_133_magnetization	magnetization	Ni2MnGa_sample_1	133	521942.7025306383	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_134_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_134_magnetization	magnetization	Ni2MnGa_sample_1	134	523469.4413496738	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_135_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_135_magnetization	magnetization	Ni2MnGa_sample_1	135	524996.1801687094	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_136_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_136_magnetization	magnetization	Ni2MnGa_sample_1	136	526523.0781473819	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_137_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_137_magnetization	magnetization	Ni2MnGa_sample_1	137	528049.8169664174	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_138_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_138_magnetization	magnetization	Ni2MnGa_sample_1	138	529576.5557854528	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_139_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_139_magnetization	magnetization	Ni2MnGa_sample_1	139	531103.374184307	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_140_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_140_magnetization	magnetization	Ni2MnGa_sample_1	140	532630.1925831609	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_141_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_141_magnetization	magnetization	Ni2MnGa_sample_1	141	534156.9314021964	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_142_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_142_magnetization	magnetization	Ni2MnGa_sample_1	142	535683.6702212319	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_143_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_143_magnetization	magnetization	Ni2MnGa_sample_1	143	537210.5681999046	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_144_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_144_magnetization	magnetization	Ni2MnGa_sample_1	144	538737.3070189401	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_145_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_145_magnetization	magnetization	Ni2MnGa_sample_1	145	540264.0458379755	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_146_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_146_magnetization	magnetization	Ni2MnGa_sample_1	146	541790.8642368296	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_147_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_147_magnetization	magnetization	Ni2MnGa_sample_1	147	543317.6826356837	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_148_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_148_magnetization	magnetization	Ni2MnGa_sample_1	148	544844.4214547192	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_149_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_149_magnetization	magnetization	Ni2MnGa_sample_1	149	546371.1602737546	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_150_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_150_magnetization	magnetization	Ni2MnGa_sample_1	150	547898.0582524271	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_151_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_151_magnetization	magnetization	Ni2MnGa_sample_1	151	549424.7970714627	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_152_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_152_magnetization	magnetization	Ni2MnGa_sample_1	152	550951.6154703167	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_153_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_153_magnetization	magnetization	Ni2MnGa_sample_1	153	552478.3542893522	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_154_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_154_magnetization	magnetization	Ni2MnGa_sample_1	154	554005.1726882062	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_155_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_155_magnetization	magnetization	Ni2MnGa_sample_1	155	555531.9115072417	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_156_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_156_magnetization	magnetization	Ni2MnGa_sample_1	156	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_157_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_157_magnetization	magnetization	Ni2MnGa_sample_1	157	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_158_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_158_magnetization	magnetization	Ni2MnGa_sample_1	158	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_159_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_159_magnetization	magnetization	Ni2MnGa_sample_1	159	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_160_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_160_magnetization	magnetization	Ni2MnGa_sample_1	160	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_161_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_161_magnetization	magnetization	Ni2MnGa_sample_1	161	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_162_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_162_magnetization	magnetization	Ni2MnGa_sample_1	162	557058.7299060959	calc_magnetization	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_111_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_111_martensite_content	martensite_content	wire_actuator_1	111	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_112_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_112_martensite_content	martensite_content	wire_actuator_1	112	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_113_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_113_martensite_content	martensite_content	wire_actuator_1	113	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_114_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_114_martensite_content	martensite_content	wire_actuator_1	114	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_115_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_115_martensite_content	martensite_content	wire_actuator_1	115	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_116_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_116_martensite_content	martensite_content	wire_actuator_1	116	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_117_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_117_martensite_content	martensite_content	wire_actuator_1	117	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_118_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_118_martensite_content	martensite_content	wire_actuator_1	118	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_119_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_119_martensite_content	martensite_content	wire_actuator_1	119	0.9567727288213013	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_120_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_120_martensite_content	martensite_content	wire_actuator_1	120	0.7191855733945405	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_121_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_121_martensite_content	martensite_content	wire_actuator_1	121	0.37903905220016804	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_122_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_122_martensite_content	martensite_content	wire_actuator_1	122	0.019369152030841108	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_123_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_123_martensite_content	martensite_content	wire_actuator_1	123	0	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_131_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_131_martensite_content	martensite_content	wire_actuator_1	131	1	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_132_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_132_martensite_content	martensite_content	wire_actuator_1	132	0.9427280128266052	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_133_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_133_martensite_content	martensite_content	wire_actuator_1	133	0.43973165987233864	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_134_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_134_martensite_content	martensite_content	wire_actuator_1	134	0.16843867087960235	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_135_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_135_martensite_content	martensite_content	wire_actuator_1	135	0.014529091286973939	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_136_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_136_martensite_content	martensite_content	wire_actuator_1	136	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_137_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_137_martensite_content	martensite_content	wire_actuator_1	137	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_138_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_138_martensite_content	martensite_content	wire_actuator_1	138	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_139_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_139_martensite_content	martensite_content	wire_actuator_1	139	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_140_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_140_martensite_content	martensite_content	wire_actuator_1	140	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_141_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_141_martensite_content	martensite_content	wire_actuator_1	141	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_142_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_142_martensite_content	martensite_content	wire_actuator_1	142	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_143_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_143_martensite_content	martensite_content	wire_actuator_1	143	0	calc_martensite_asc	1
cp_mod_maximal_blocking_stress_electrostatic_pressure_model_elastomer_1_0_maximal_blocking_stress	mod_electrostatic_pressure_model_maximal_blocking_stress_elastomer_1_0_maximal_blocking_stress	maximal_blocking_stress	elastomer_1	0	277689.7728	electrostatic_pressure_model	1
cp_mod_maximum_strain_mag_calc_maximum_strain_mag_Ni2MnGa_sample_1_0_maximum_strain_mag	mod_calc_maximum_strain_mag_maximum_strain_mag_Ni2MnGa_sample_1_0_maximum_strain_mag	maximum_strain_mag	Ni2MnGa_sample_1	0	0.06779661016949146	calc_maximum_strain_mag	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_110_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_110_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	110	0	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_111_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_111_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	111	267550000	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_112_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_112_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	112	471068000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_113_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_113_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	113	471818000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_114_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_114_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	114	472568000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_115_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_115_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	115	473318000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_116_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_116_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	116	474068000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_117_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_117_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	117	474818000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_118_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_118_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	118	475568000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_119_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_119_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	119	476318000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_120_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_120_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	120	477068000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_121_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_121_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	121	477818000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v3_wire_actuator_NiTi#6_122_mechanic_stress	mod_1d_non_linear_mech_mod_v3_mechanic_stress_wire_actuator_NiTi#6_122_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	122	525745000	1d_non_linear_mech_mod_v3	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_210_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_210_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	210	0	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_211_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_211_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	211	190775000	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_212_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_212_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	212	335944000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_213_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_213_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	213	336694000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_214_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_214_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	214	337444000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_215_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_215_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	215	338194000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_216_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_216_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	216	338944000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_217_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_217_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	217	339694000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_218_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_218_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	218	340444000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_219_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_219_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	219	341194000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_220_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_220_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	220	341944000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_221_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_221_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	221	342694000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v3_wire_actuator_NiTi#6_222_mechanic_stress	mod_1d_non_linear_mech_mod_v3_mechanic_stress_wire_actuator_NiTi#6_222_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	222	390621000	1d_non_linear_mech_mod_v3	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_310_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_310_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	310	0	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_311_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_311_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	311	114000000	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_312_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_312_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	312	200820000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_313_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_313_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	313	201570000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_314_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_314_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	314	202320000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_315_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_315_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	315	203070000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_316_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_316_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	316	203820000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_317_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_317_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	317	204570000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_318_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_318_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	318	205320000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_319_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_319_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	319	206070000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_320_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_320_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	320	206820000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_321_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_321_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	321	207570000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v3_wire_actuator_NiTi#6_322_mechanic_stress	mod_1d_non_linear_mech_mod_v3_mechanic_stress_wire_actuator_NiTi#6_322_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	322	255496999.99999997	1d_non_linear_mech_mod_v3	1
cp_mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_1	mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_1	Yeoh_vector_1	sample1_elastomer1	0	194409.9644002127	calc_yeoh_minR	1
cp_mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_2	mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_2	Yeoh_vector_2	sample1_elastomer1	0	-2275.4747645269717	calc_yeoh_minR	1
cp_mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_3	mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_3	Yeoh_vector_3	sample1_elastomer1	0	527.8018381989875	calc_yeoh_minR	1
cp_mod_calc_young_modulus_minR_sample1_elastomer1young_modulus	mod_calc_young_modulus_minR_sample1_elastomer1young_modulus	young_modulus	sample1_elastomer1	0	652412.0602756739	calc_young_modulus_minR	1
cp_mod_calc_young_modulus_Neo_Hookean_minR_sample1_elastomer1young_modulus_neo_Hookean	mod_calc_young_modulus_Neo_Hookean_minR_sample1_elastomer1young_modulus_neo_Hookean	young_modulus_neo_Hookean	sample1_elastomer1	0	1178216.3521274924	calc_young_modulus_Neo_Hookean_minR	1
cp_mod_A00_matrix_model_1_M420_0_1	mod_matrix_model_1_A00_M420_0	A00	M420	0	1.54e-11	matrix_model_1	2
cp_mod_A00_matrix_model_1_M420_0_2	mod_matrix_model_1_A00_M420_0	A00	M420	0	-5.70000000000003e-12	matrix_model_1	2
cp_mod_A00_matrix_model_1_M420_0_3	mod_matrix_model_1_A00_M420_0	A00	M420	0	-6.49999999999999e-12	matrix_model_1	2
cp_mod_A00_matrix_model_1_M420_0_4	mod_matrix_model_1_A00_M420_0	A00	M420	0	1.87e-11	matrix_model_1	2
cp_mod_A00_matrix_model_1_M420_0_5	mod_matrix_model_1_A00_M420_0	A00	M420	0	4.5e-11	matrix_model_1	2
cp_mod_A11_matrix_model_1_M420_0_1	mod_matrix_model_1_A11_M420_0	A11	M420	0	127420865309.707	matrix_model_1	2
cp_mod_A11_matrix_model_1_M420_0_2	mod_matrix_model_1_A11_M420_0	A11	M420	0	80027500380.7973	matrix_model_1	2
cp_mod_A11_matrix_model_1_M420_0_3	mod_matrix_model_1_A11_M420_0	A11	M420	0	52698431980.2742	matrix_model_1	2
cp_mod_A11_matrix_model_1_M420_0_4	mod_matrix_model_1_A11_M420_0	A11	M420	0	128773831986.699	matrix_model_1	2
cp_mod_A11_matrix_model_1_M420_0_5	mod_matrix_model_1_A11_M420_0	A11	M420	0	39148517532.8784	matrix_model_1	2
cp_mod_B00_matrix_model_2_M420_0_3	mod_matrix_model_2_B00_M420_0	B00	M420	0	5.25000000000002e-10	matrix_model_2	2
cp_mod_B00_matrix_model_2_M420_0_1	mod_matrix_model_2_B00_M420_0	B00	M420	0	-1.6e-10	matrix_model_2	2
cp_mod_B00_matrix_model_2_M420_0_2	mod_matrix_model_2_B00_M420_0	B00	M420	0	3.55e-10	matrix_model_2	2
cp_mod_B11_matrix_model_2_M420_0_3	mod_matrix_model_2_B11_M420_0	B11	M420	0	1450825312.34196	matrix_model_2	2
cp_mod_B11_matrix_model_2_M420_0_1	mod_matrix_model_2_B11_M420_0	B11	M420	0	-964471223.088693	matrix_model_2	2
cp_mod_B11_matrix_model_2_M420_0_2	mod_matrix_model_2_B11_M420_0	B11	M420	0	2059980296.10146	matrix_model_2	2
cp_mod_C00_matrix_model_3_M420_0_1	mod_matrix_model_3_C00_M420_0	C00	M420	0	1.41664e-08	matrix_model_3	2
cp_mod_C00_matrix_model_3_M420_0_2	mod_matrix_model_3_C00_M420_0	C00	M420	0	1.41664e-08	matrix_model_3	2
cp_mod_C11_matrix_model_3_M420_0_1	mod_matrix_model_3_C11_M420_0	C11	M420	0	124356455.343597	matrix_model_3	2
cp_mod_C11_matrix_model_3_M420_0_2	mod_matrix_model_3_C11_M420_0	C11	M420	0	144719941.118566	matrix_model_3	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_1_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_1_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	1	0	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_2_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_2_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	2	0.10026755000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_3_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_3_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	3	0.12479048999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_4_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_4_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	4	0.1356435	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_5_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_5_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	5	0.1464965	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_6_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_6_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	6	0.15734951	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_7_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_7_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	7	0.16820252000000002	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_8_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_8_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	8	0.17905554	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_9_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_9_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	9	0.18990854000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_10_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_10_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	10	0.20076155999999998	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_11_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_11_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	11	0.21161454999999998	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_12_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_12_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	12	0.22246757	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_13_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_13_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	13	0.23332060000000002	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_14_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_14_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	14	0.24417360000000002	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_15_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_15_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	15	0.2550266	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_16_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_16_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	16	0.2658796	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_17_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_17_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	17	0.27673263000000003	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_18_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_18_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	18	0.28758564999999997	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_19_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_19_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	19	0.29843860000000005	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_20_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_20_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	20	0.30929164	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_21_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_21_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	21	0.32014465999999997	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_22_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_22_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	22	0.3309977	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_23_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_23_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	23	0.3418507	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_24_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_24_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	24	0.3527036700000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_25_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_25_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	25	0.36355670000000007	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_26_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_26_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	26	0.37440970000000007	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_27_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_27_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	27	0.38526273000000005	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_28_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_28_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	28	0.3961157600000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_29_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_29_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	29	0.40696870000000007	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_30_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_30_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	30	0.41782173999999994	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_31_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_31_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	31	0.42867476000000004	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_32_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_32_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	32	0.43952779999999997	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_33_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_33_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	33	0.45038075	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_34_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_34_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	34	0.46123379999999997	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_35_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_35_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	35	0.4720868	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_36_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_36_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	36	0.4829398	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_37_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_37_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	37	0.49379284	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_38_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_38_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	38	0.5046457999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_39_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_39_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	39	0.5154987999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_40_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_40_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	40	0.5263517999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_41_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_41_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	41	0.5372049000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_42_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_42_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	42	0.5480578999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_43_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_43_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	43	0.5589109	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_44_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_44_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	44	0.5697639	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_45_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_45_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	45	0.5806168	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_46_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_46_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	46	0.5914699	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_47_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_47_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	47	0.6023229000000002	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_48_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_48_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	48	0.6131759	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_49_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_49_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	49	0.6240289	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_50_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_50_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	50	0.634882	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_51_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_51_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	51	0.6457350000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_52_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_52_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	52	0.6565880000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_53_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_53_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	53	0.667441	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_54_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_54_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	54	0.6782939000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_55_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_55_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	55	0.689147	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_56_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_56_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	56	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_57_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_57_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	57	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_58_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_58_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	58	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_59_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_59_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	59	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_60_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_60_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	60	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_61_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_61_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	61	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_62_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_62_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	62	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_101_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_101_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	101	0	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_102_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_102_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	102	0.4230082	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_103_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_103_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	103	0.5983164300000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_104_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_104_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	104	0.60023496	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_105_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_105_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	105	0.60215354	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_106_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_106_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	106	0.6040721	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_107_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_107_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	107	0.60599065	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_108_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_108_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	108	0.60790923	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_109_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_109_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	109	0.6098277599999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_110_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_110_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	110	0.6117463399999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_111_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_111_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	111	0.61366487	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_112_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_112_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	112	0.61558345	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_113_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_113_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	113	0.617502	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_114_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_114_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	114	0.61942056	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_115_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_115_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	115	0.62133915	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_116_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_116_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	116	0.6232576700000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_117_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_117_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	117	0.6251762500000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_118_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_118_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	118	0.6270948	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_119_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_119_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	119	0.62901336	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_120_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_120_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	120	0.6309319000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_121_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_121_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	121	0.6328505	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_122_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_122_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	122	0.63476906	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_123_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_123_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	123	0.6366876	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_124_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_124_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	124	0.63860617	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_125_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_125_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	125	0.6405247000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_126_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_126_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	126	0.6424433	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_127_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_127_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	127	0.6443618	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_128_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_128_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	128	0.6462803999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_129_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_129_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	129	0.648199	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_130_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_130_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	130	0.6501175	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_131_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_131_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	131	0.652036	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_132_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_132_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	132	0.6539547000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_133_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_133_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	133	0.6558732	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_134_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_134_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	134	0.6577917000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_135_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_135_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	135	0.6597102000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_136_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_136_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	136	0.6616289000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_137_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_137_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	137	0.6635474	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_138_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_138_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	138	0.6654659	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_139_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_139_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	139	0.6673845	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_140_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_140_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	140	0.6693031	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_141_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_141_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	141	0.6712216	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_142_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_142_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	142	0.6731400999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_143_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_143_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	143	0.6750588000000002	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_144_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_144_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	144	0.6769773000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_145_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_145_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	145	0.6788957999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_146_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_146_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	146	0.6808144	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_147_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_147_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	147	0.6827330000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_148_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_148_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	148	0.6846515000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_149_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_149_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	149	0.68657	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_150_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_150_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	150	0.6884887	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_151_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_151_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	151	0.6904072	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_152_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_152_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	152	0.6923258	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_153_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_153_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	153	0.6942442999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_154_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_154_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	154	0.6961628999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_155_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_155_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	155	0.6980813999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_156_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_156_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	156	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_157_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_157_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	157	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_158_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_158_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	158	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_159_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_159_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	159	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_160_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_160_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	160	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_161_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_161_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	161	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_162_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_162_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	162	0.7000000000000001	calc_mag_polarization	2
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_1	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	1	0	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_2	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	2	501.33775	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_3	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	3	1626.62795	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_4	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	4	2928.7979	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_5	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	5	4339.4979	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_6	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	6	5858.72795	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_7	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	7	7486.4881000000005	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_8	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	8	9222.778400000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_9	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	9	11067.598800000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_10	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	10	13020.949300000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_11	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	11	15082.829850000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_12	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	12	17253.24045	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_13	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	13	19532.1813	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_14	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	14	21919.6523	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_15	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	15	24415.6533	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_16	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	16	27020.1843	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_17	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	17	29733.245450000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_18	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	18	32554.836850000003	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_19	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	19	35484.9581	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_20	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	20	38523.609300000004	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_21	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	21	41670.7908	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_22	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	22	44926.5026	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_23	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	23	48290.7446	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_24	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	24	51763.516449999996	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_25	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	25	55344.8183	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_26	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	26	59034.6503	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_27	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	27	62833.01245	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_28	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	28	66739.90490000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_29	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	29	70755.32720000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_30	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	30	74879.27940000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_31	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	31	79111.76190000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_32	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	32	83452.77470000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_33	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	33	87902.31745	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_34	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	34	92460.3902	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_35	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	35	97126.9932	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_36	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	36	101902.1262	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_37	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	37	106785.7894	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_38	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	38	111777.98259999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_39	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	39	116878.70559999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_40	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	40	122087.95859999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_41	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	41	127405.74209999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_42	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	42	132832.0561	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_43	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	43	138366.9001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_44	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	44	144010.2741	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_45	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	45	149762.1776	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_46	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	46	155622.61109999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_47	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	47	161591.5751	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_48	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	48	167669.0691	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_49	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	49	173855.0931	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_50	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	50	180149.6476	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_51	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	51	186552.7326	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_52	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	52	193064.34759999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_53	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	53	199684.49259999997	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_54	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	54	206413.16709999996	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_55	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	55	213250.37159999995	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_56	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	56	220196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_57	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	57	227196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_58	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	58	234196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_59	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	59	241196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_60	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	60	248196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_61	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	61	255196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_62	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	62	262196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_101	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	101	0	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_102	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	102	2115.041	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_103	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	103	7221.6641500000005	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_104	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	104	13214.4211	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_105	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	105	19226.3636	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_106	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	106	25257.4918	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_107	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	107	31307.80555	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_108	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	108	37377.304950000005	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_109	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	109	43465.98990000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_110	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	110	49573.860400000005	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_111	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	111	55700.916450000004	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_112	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	112	61847.158050000005	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_113	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	113	68012.5853	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_114	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	114	74197.19810000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_115	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	115	80400.99665000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_116	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	116	86623.98075000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_117	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	117	92866.15035000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_118	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	118	99127.5056	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_119	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	119	105408.0464	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_120	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	120	111707.7727	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_121	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	121	118026.6847	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_122	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	122	124364.7825	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_123	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	123	130722.0658	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_124	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	124	137098.53465	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_125	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	125	143494.18899999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_126	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	126	149909.02899999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_127	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	127	156343.05449999997	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_128	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	128	162796.26549999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_129	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	129	169268.66249999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_130	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	130	175760.24499999997	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_131	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	131	182271.01249999995	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_132	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	132	188800.96599999996	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_133	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	133	195350.10549999995	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_134	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	134	201918.42999999993	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_135	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	135	208505.93949999995	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_136	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	136	215112.63499999995	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_137	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	137	221738.51649999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_138	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	138	228383.58299999993	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_139	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	139	235047.83499999993	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_140	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	140	241731.27299999993	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_141	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	141	248433.89649999992	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_142	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	142	255155.7049999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_143	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	143	261896.6994999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_144	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	144	268656.8799999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_145	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	145	275436.2454999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_146	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	146	282234.7964999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_147	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	147	289052.5334999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_148	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	148	295889.4559999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_149	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	149	302745.5634999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_150	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	150	309620.85699999984	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_151	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	151	316515.33649999986	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_152	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	152	323429.00149999984	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_153	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	153	330361.85199999984	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_154	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	154	337313.88799999986	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_155	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	155	344285.10949999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_156	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	156	351275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_157	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	157	358275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_158	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	158	365275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_159	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	159	372275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_160	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	160	379275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_161	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	161	386275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_162	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	162	393275.51649999985	calc_mag_energy_density	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_2_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_2_relative_permeability	relative_permeability	Ni2MnGa_sample_1	2	8.979273436256566	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_3_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_3_relative_permeability	relative_permeability	Ni2MnGa_sample_1	3	5.965402275982811	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_4_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_4_relative_permeability	relative_permeability	Ni2MnGa_sample_1	4	4.598161706191309	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_5_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_5_relative_permeability	relative_permeability	Ni2MnGa_sample_1	5	3.9145412223460134	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_6_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_6_relative_permeability	relative_permeability	Ni2MnGa_sample_1	6	3.504369091198472	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_7_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_7_relative_permeability	relative_permeability	Ni2MnGa_sample_1	7	3.2309210037667784	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_8_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_8_relative_permeability	relative_permeability	Ni2MnGa_sample_1	8	3.0356010550010235	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_9_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_9_relative_permeability	relative_permeability	Ni2MnGa_sample_1	9	2.889110894477161	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_10_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_10_relative_permeability	relative_permeability	Ni2MnGa_sample_1	10	2.7751742798026418	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_11_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_11_relative_permeability	relative_permeability	Ni2MnGa_sample_1	11	2.6840247493235716	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_12_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_12_relative_permeability	relative_permeability	Ni2MnGa_sample_1	12	2.6094480777856557	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_13_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_13_relative_permeability	relative_permeability	Ni2MnGa_sample_1	13	2.547300917820574	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_14_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_14_relative_permeability	relative_permeability	Ni2MnGa_sample_1	14	2.494714675742847	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_15_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_15_relative_permeability	relative_permeability	Ni2MnGa_sample_1	15	2.449640753961938	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_16_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_16_relative_permeability	relative_permeability	Ni2MnGa_sample_1	16	2.410576688418484	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_17_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_17_relative_permeability	relative_permeability	Ni2MnGa_sample_1	17	2.376395780280121	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_18_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_18_relative_permeability	relative_permeability	Ni2MnGa_sample_1	18	2.346236108640496	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_19_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_19_relative_permeability	relative_permeability	Ni2MnGa_sample_1	19	2.3194272021504236	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_20_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_20_relative_permeability	relative_permeability	Ni2MnGa_sample_1	20	2.295440662774236	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_21_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_21_relative_permeability	relative_permeability	Ni2MnGa_sample_1	21	2.2738526977558493	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_22_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_22_relative_permeability	relative_permeability	Ni2MnGa_sample_1	22	2.2543208051961834	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_23_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_23_relative_permeability	relative_permeability	Ni2MnGa_sample_1	23	2.2365643945422713	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_24_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_24_relative_permeability	relative_permeability	Ni2MnGa_sample_1	24	2.2203519157976324	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_25_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_25_relative_permeability	relative_permeability	Ni2MnGa_sample_1	25	2.205490675897926	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_26_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_26_relative_permeability	relative_permeability	Ni2MnGa_sample_1	26	2.1918182396944137	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_27_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_27_relative_permeability	relative_permeability	Ni2MnGa_sample_1	27	2.179197621175578	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_28_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_28_relative_permeability	relative_permeability	Ni2MnGa_sample_1	28	2.167511863287767	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_29_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_29_relative_permeability	relative_permeability	Ni2MnGa_sample_1	29	2.156660546599668	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_30_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_30_relative_permeability	relative_permeability	Ni2MnGa_sample_1	30	2.1465578709928814	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_31_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_31_relative_permeability	relative_permeability	Ni2MnGa_sample_1	31	2.1371286540400023	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_32_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_32_relative_permeability	relative_permeability	Ni2MnGa_sample_1	32	2.128307825006546	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_33_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_33_relative_permeability	relative_permeability	Ni2MnGa_sample_1	33	2.1200380739694413	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_34_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_34_relative_permeability	relative_permeability	Ni2MnGa_sample_1	34	2.1122697611158534	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_35_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_35_relative_permeability	relative_permeability	Ni2MnGa_sample_1	35	2.104958290812744	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_36_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_36_relative_permeability	relative_permeability	Ni2MnGa_sample_1	36	2.098064618812669	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_37_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_37_relative_permeability	relative_permeability	Ni2MnGa_sample_1	37	2.0915540170123967	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_38_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_38_relative_permeability	relative_permeability	Ni2MnGa_sample_1	38	2.0853951675692883	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_39_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_39_relative_permeability	relative_permeability	Ni2MnGa_sample_1	39	2.0795605518650992	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_40_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_40_relative_permeability	relative_permeability	Ni2MnGa_sample_1	40	2.074025147222664	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_41_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_41_relative_permeability	relative_permeability	Ni2MnGa_sample_1	41	2.0687667117618975	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_42_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_42_relative_permeability	relative_permeability	Ni2MnGa_sample_1	42	2.0637645912508784	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_43_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_43_relative_permeability	relative_permeability	Ni2MnGa_sample_1	43	2.0590006669546694	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_44_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_44_relative_permeability	relative_permeability	Ni2MnGa_sample_1	44	2.054458320532704	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_45_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_45_relative_permeability	relative_permeability	Ni2MnGa_sample_1	45	2.0501222635394214	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_46_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_46_relative_permeability	relative_permeability	Ni2MnGa_sample_1	46	2.0459792738783666	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_47_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_47_relative_permeability	relative_permeability	Ni2MnGa_sample_1	47	2.0420162412029708	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_48_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_48_relative_permeability	relative_permeability	Ni2MnGa_sample_1	48	2.0382218482158883	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_49_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_49_relative_permeability	relative_permeability	Ni2MnGa_sample_1	49	2.0345855549366014	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_50_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_50_relative_permeability	relative_permeability	Ni2MnGa_sample_1	50	2.031097844198956	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_51_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_51_relative_permeability	relative_permeability	Ni2MnGa_sample_1	51	2.0277494827311795	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_52_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_52_relative_permeability	relative_permeability	Ni2MnGa_sample_1	52	2.024532429556257	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_53_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_53_relative_permeability	relative_permeability	Ni2MnGa_sample_1	53	2.0214391091957538	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_54_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_54_relative_permeability	relative_permeability	Ni2MnGa_sample_1	54	2.018462367754858	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_55_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_55_relative_permeability	relative_permeability	Ni2MnGa_sample_1	55	2.015596170737027	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_56_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_56_relative_permeability	relative_permeability	Ni2MnGa_sample_1	56	2.0128340543747196	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_57_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_57_relative_permeability	relative_permeability	Ni2MnGa_sample_1	57	1.9947477319751712	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_58_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_58_relative_permeability	relative_permeability	Ni2MnGa_sample_1	58	1.9772960173791156	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_59_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_59_relative_permeability	relative_permeability	Ni2MnGa_sample_1	59	1.9604460860449928	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_60_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_60_relative_permeability	relative_permeability	Ni2MnGa_sample_1	60	1.9441673388238914	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_61_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_61_relative_permeability	relative_permeability	Ni2MnGa_sample_1	61	1.9284312165101598	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_62_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_62_relative_permeability	relative_permeability	Ni2MnGa_sample_1	62	1.9132110326329441	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_102_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_102_relative_permeability	relative_permeability	Ni2MnGa_sample_1	102	34.66291580455197	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_103_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_103_relative_permeability	relative_permeability	Ni2MnGa_sample_1	103	24.806956469839257	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_104_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_104_relative_permeability	relative_permeability	Ni2MnGa_sample_1	104	16.9221964029922	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_105_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_105_relative_permeability	relative_permeability	Ni2MnGa_sample_1	105	12.97981736431641	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_106_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_106_relative_permeability	relative_permeability	Ni2MnGa_sample_1	106	10.614389622791661	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_107_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_107_relative_permeability	relative_permeability	Ni2MnGa_sample_1	107	9.037437662475462	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_108_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_108_relative_permeability	relative_permeability	Ni2MnGa_sample_1	108	7.911043746163116	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_109_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_109_relative_permeability	relative_permeability	Ni2MnGa_sample_1	109	7.066247811554989	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_110_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_110_relative_permeability	relative_permeability	Ni2MnGa_sample_1	110	6.409184748969883	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_111_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_111_relative_permeability	relative_permeability	Ni2MnGa_sample_1	111	5.883533901002706	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_112_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_112_relative_permeability	relative_permeability	Ni2MnGa_sample_1	112	5.453456296210554	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_113_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_113_relative_permeability	relative_permeability	Ni2MnGa_sample_1	113	5.095058093267547	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_114_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_114_relative_permeability	relative_permeability	Ni2MnGa_sample_1	114	4.791798136607941	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_115_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_115_relative_permeability	relative_permeability	Ni2MnGa_sample_1	115	4.531861201427889	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_116_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_116_relative_permeability	relative_permeability	Ni2MnGa_sample_1	116	4.306582152899359	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_117_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_117_relative_permeability	relative_permeability	Ni2MnGa_sample_1	117	4.1094632838612135	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_118_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_118_relative_permeability	relative_permeability	Ni2MnGa_sample_1	118	3.9355347295690515	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_119_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_119_relative_permeability	relative_permeability	Ni2MnGa_sample_1	119	3.780931614409252	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_120_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_120_relative_permeability	relative_permeability	Ni2MnGa_sample_1	120	3.6426024276033075	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_121_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_121_relative_permeability	relative_permeability	Ni2MnGa_sample_1	121	3.5181063982174123	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_122_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_122_relative_permeability	relative_permeability	Ni2MnGa_sample_1	122	3.4054669819543286	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_123_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_123_relative_permeability	relative_permeability	Ni2MnGa_sample_1	123	3.303067440278964	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_124_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_124_relative_permeability	relative_permeability	Ni2MnGa_sample_1	124	3.2095723103751324	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_125_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_125_relative_permeability	relative_permeability	Ni2MnGa_sample_1	125	3.1238683086635897	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_126_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_126_relative_permeability	relative_permeability	Ni2MnGa_sample_1	126	3.0450208499124622	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_127_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_127_relative_permeability	relative_permeability	Ni2MnGa_sample_1	127	2.9722382742198117	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_128_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_128_relative_permeability	relative_permeability	Ni2MnGa_sample_1	128	2.9048472951703888	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_129_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_129_relative_permeability	relative_permeability	Ni2MnGa_sample_1	129	2.8422699574816397	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_130_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_130_relative_permeability	relative_permeability	Ni2MnGa_sample_1	130	2.784008023841016	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_131_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_131_relative_permeability	relative_permeability	Ni2MnGa_sample_1	131	2.729630219109767	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_132_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_132_relative_permeability	relative_permeability	Ni2MnGa_sample_1	132	2.6787611732632346	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_133_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_133_relative_permeability	relative_permeability	Ni2MnGa_sample_1	133	2.6310709454082444	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_134_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_134_relative_permeability	relative_permeability	Ni2MnGa_sample_1	134	2.5862710343929507	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_135_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_135_relative_permeability	relative_permeability	Ni2MnGa_sample_1	135	2.54410641226091	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_136_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_136_relative_permeability	relative_permeability	Ni2MnGa_sample_1	136	2.504351651849663	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_137_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_137_relative_permeability	relative_permeability	Ni2MnGa_sample_1	137	2.466805047128937	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_138_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_138_relative_permeability	relative_permeability	Ni2MnGa_sample_1	138	2.4312879886093324	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_139_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_139_relative_permeability	relative_permeability	Ni2MnGa_sample_1	139	2.397640458379755	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_140_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_140_relative_permeability	relative_permeability	Ni2MnGa_sample_1	140	2.3657184425209254	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_141_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_141_relative_permeability	relative_permeability	Ni2MnGa_sample_1	141	2.335392328505491	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_142_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_142_relative_permeability	relative_permeability	Ni2MnGa_sample_1	142	2.306545537124956	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_143_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_143_relative_permeability	relative_permeability	Ni2MnGa_sample_1	143	2.2790727814283445	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_144_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_144_relative_permeability	relative_permeability	Ni2MnGa_sample_1	144	2.2528774581835815	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_145_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_145_relative_permeability	relative_permeability	Ni2MnGa_sample_1	145	2.227872831449944	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_146_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_146_relative_permeability	relative_permeability	Ni2MnGa_sample_1	146	2.2039796983040656	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_147_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_147_relative_permeability	relative_permeability	Ni2MnGa_sample_1	147	2.181125397034095	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_148_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_148_relative_permeability	relative_permeability	Ni2MnGa_sample_1	148	2.159243449903658	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_149_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_149_relative_permeability	relative_permeability	Ni2MnGa_sample_1	149	2.1382732505703217	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_150_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_150_relative_permeability	relative_permeability	Ni2MnGa_sample_1	150	2.1181593025559735	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_151_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_151_relative_permeability	relative_permeability	Ni2MnGa_sample_1	151	2.0988495941429255	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_152_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_152_relative_permeability	relative_permeability	Ni2MnGa_sample_1	152	2.080297285235915	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_153_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_153_relative_permeability	relative_permeability	Ni2MnGa_sample_1	153	2.0624583736333695	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_154_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_154_relative_permeability	relative_permeability	Ni2MnGa_sample_1	154	2.045292778656993	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_155_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_155_relative_permeability	relative_permeability	Ni2MnGa_sample_1	155	2.0287627990874846	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_156_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_156_relative_permeability	relative_permeability	Ni2MnGa_sample_1	156	2.0128340543747196	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_157_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_157_relative_permeability	relative_permeability	Ni2MnGa_sample_1	157	1.9947477319751712	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_158_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_158_relative_permeability	relative_permeability	Ni2MnGa_sample_1	158	1.9772960173791156	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_159_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_159_relative_permeability	relative_permeability	Ni2MnGa_sample_1	159	1.9604460860449928	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_160_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_160_relative_permeability	relative_permeability	Ni2MnGa_sample_1	160	1.9441673388238914	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_161_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_161_relative_permeability	relative_permeability	Ni2MnGa_sample_1	161	1.9284312165101598	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_162_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_162_relative_permeability	relative_permeability	Ni2MnGa_sample_1	162	1.9132110326329441	calc_rel_perm	3
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_2_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_2_flux_density	flux_density	Ni2MnGa_sample_1	2	0.11283355000000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_3_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_3_flux_density	flux_density	Ni2MnGa_sample_1	3	0.14992249	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_4_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_4_flux_density	flux_density	Ni2MnGa_sample_1	4	0.17334149999999995	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_5_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_5_flux_density	flux_density	Ni2MnGa_sample_1	5	0.19676050000000003	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_6_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_6_flux_density	flux_density	Ni2MnGa_sample_1	6	0.22017951	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_7_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_7_flux_density	flux_density	Ni2MnGa_sample_1	7	0.24359851999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_8_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_8_flux_density	flux_density	Ni2MnGa_sample_1	8	0.26701754	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_9_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_9_flux_density	flux_density	Ni2MnGa_sample_1	9	0.29043654	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_10_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_10_flux_density	flux_density	Ni2MnGa_sample_1	10	0.31385556	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_11_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_11_flux_density	flux_density	Ni2MnGa_sample_1	11	0.33727455	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_12_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_12_flux_density	flux_density	Ni2MnGa_sample_1	12	0.3606935700000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_13_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_13_flux_density	flux_density	Ni2MnGa_sample_1	13	0.3841126	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_14_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_14_flux_density	flux_density	Ni2MnGa_sample_1	14	0.4075316	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_15_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_15_flux_density	flux_density	Ni2MnGa_sample_1	15	0.4309505999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_16_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_16_flux_density	flux_density	Ni2MnGa_sample_1	16	0.4543696	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_17_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_17_flux_density	flux_density	Ni2MnGa_sample_1	17	0.47778862999999994	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_18_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_18_flux_density	flux_density	Ni2MnGa_sample_1	18	0.50120765	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_19_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_19_flux_density	flux_density	Ni2MnGa_sample_1	19	0.5246266	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_20_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_20_flux_density	flux_density	Ni2MnGa_sample_1	20	0.54804564	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_21_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_21_flux_density	flux_density	Ni2MnGa_sample_1	21	0.57146466	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_22_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_22_flux_density	flux_density	Ni2MnGa_sample_1	22	0.5948837	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_23_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_23_flux_density	flux_density	Ni2MnGa_sample_1	23	0.6183027	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_24_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_24_flux_density	flux_density	Ni2MnGa_sample_1	24	0.6417216700000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_25_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_25_flux_density	flux_density	Ni2MnGa_sample_1	25	0.6651407	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_26_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_26_flux_density	flux_density	Ni2MnGa_sample_1	26	0.6885597	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_27_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_27_flux_density	flux_density	Ni2MnGa_sample_1	27	0.7119787300000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_28_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_28_flux_density	flux_density	Ni2MnGa_sample_1	28	0.7353977600000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_29_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_29_flux_density	flux_density	Ni2MnGa_sample_1	29	0.7588167000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_30_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_30_flux_density	flux_density	Ni2MnGa_sample_1	30	0.78223574	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_31_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_31_flux_density	flux_density	Ni2MnGa_sample_1	31	0.80565476	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_32_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_32_flux_density	flux_density	Ni2MnGa_sample_1	32	0.8290738	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_33_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_33_flux_density	flux_density	Ni2MnGa_sample_1	33	0.8524927499999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_34_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_34_flux_density	flux_density	Ni2MnGa_sample_1	34	0.8759117999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_35_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_35_flux_density	flux_density	Ni2MnGa_sample_1	35	0.8993308	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_36_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_36_flux_density	flux_density	Ni2MnGa_sample_1	36	0.9227498	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_37_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_37_flux_density	flux_density	Ni2MnGa_sample_1	37	0.94616884	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_38_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_38_flux_density	flux_density	Ni2MnGa_sample_1	38	0.9695878000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_39_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_39_flux_density	flux_density	Ni2MnGa_sample_1	39	0.9930067999999997	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_40_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_40_flux_density	flux_density	Ni2MnGa_sample_1	40	1.0164257999999997	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_41_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_41_flux_density	flux_density	Ni2MnGa_sample_1	41	1.0398449	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_42_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_42_flux_density	flux_density	Ni2MnGa_sample_1	42	1.0632639000000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_43_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_43_flux_density	flux_density	Ni2MnGa_sample_1	43	1.0866828999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_44_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_44_flux_density	flux_density	Ni2MnGa_sample_1	44	1.1101019	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_45_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_45_flux_density	flux_density	Ni2MnGa_sample_1	45	1.1335208	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_46_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_46_flux_density	flux_density	Ni2MnGa_sample_1	46	1.1569399	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_47_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_47_flux_density	flux_density	Ni2MnGa_sample_1	47	1.1803589000000003	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_48_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_48_flux_density	flux_density	Ni2MnGa_sample_1	48	1.2037779000000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_49_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_49_flux_density	flux_density	Ni2MnGa_sample_1	49	1.2271969	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_50_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_50_flux_density	flux_density	Ni2MnGa_sample_1	50	1.250616	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_51_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_51_flux_density	flux_density	Ni2MnGa_sample_1	51	1.274035	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_52_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_52_flux_density	flux_density	Ni2MnGa_sample_1	52	1.297454	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_53_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_53_flux_density	flux_density	Ni2MnGa_sample_1	53	1.3208729999999997	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_54_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_54_flux_density	flux_density	Ni2MnGa_sample_1	54	1.3442919	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_55_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_55_flux_density	flux_density	Ni2MnGa_sample_1	55	1.3677110000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_56_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_56_flux_density	flux_density	Ni2MnGa_sample_1	56	1.39113	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_57_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_57_flux_density	flux_density	Ni2MnGa_sample_1	57	1.4036959999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_58_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_58_flux_density	flux_density	Ni2MnGa_sample_1	58	1.4162620000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_59_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_59_flux_density	flux_density	Ni2MnGa_sample_1	59	1.428828	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_60_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_60_flux_density	flux_density	Ni2MnGa_sample_1	60	1.441394	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_61_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_61_flux_density	flux_density	Ni2MnGa_sample_1	61	1.4539600000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_62_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_62_flux_density	flux_density	Ni2MnGa_sample_1	62	1.466526	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_102_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_102_flux_density	flux_density	Ni2MnGa_sample_1	102	0.4355742000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_103_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_103_flux_density	flux_density	Ni2MnGa_sample_1	103	0.6234484300000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_104_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_104_flux_density	flux_density	Ni2MnGa_sample_1	104	0.6379329599999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_105_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_105_flux_density	flux_density	Ni2MnGa_sample_1	105	0.6524175400000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_106_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_106_flux_density	flux_density	Ni2MnGa_sample_1	106	0.6669021	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_107_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_107_flux_density	flux_density	Ni2MnGa_sample_1	107	0.68138665	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_108_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_108_flux_density	flux_density	Ni2MnGa_sample_1	108	0.69587123	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_109_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_109_flux_density	flux_density	Ni2MnGa_sample_1	109	0.7103557599999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_110_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_110_flux_density	flux_density	Ni2MnGa_sample_1	110	0.7248403399999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_111_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_111_flux_density	flux_density	Ni2MnGa_sample_1	111	0.73932487	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_112_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_112_flux_density	flux_density	Ni2MnGa_sample_1	112	0.7538094500000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_113_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_113_flux_density	flux_density	Ni2MnGa_sample_1	113	0.768294	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_114_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_114_flux_density	flux_density	Ni2MnGa_sample_1	114	0.7827785599999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_115_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_115_flux_density	flux_density	Ni2MnGa_sample_1	115	0.7972631499999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_116_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_116_flux_density	flux_density	Ni2MnGa_sample_1	116	0.81174767	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_117_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_117_flux_density	flux_density	Ni2MnGa_sample_1	117	0.8262322500000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_118_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_118_flux_density	flux_density	Ni2MnGa_sample_1	118	0.8407167999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_119_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_119_flux_density	flux_density	Ni2MnGa_sample_1	119	0.8552013599999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_120_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_120_flux_density	flux_density	Ni2MnGa_sample_1	120	0.8696859	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_121_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_121_flux_density	flux_density	Ni2MnGa_sample_1	121	0.8841705000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_122_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_122_flux_density	flux_density	Ni2MnGa_sample_1	122	0.8986550600000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_123_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_123_flux_density	flux_density	Ni2MnGa_sample_1	123	0.9131396	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_124_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_124_flux_density	flux_density	Ni2MnGa_sample_1	124	0.92762417	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_125_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_125_flux_density	flux_density	Ni2MnGa_sample_1	125	0.9421086999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_126_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_126_flux_density	flux_density	Ni2MnGa_sample_1	126	0.9565932999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_127_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_127_flux_density	flux_density	Ni2MnGa_sample_1	127	0.9710778	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_128_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_128_flux_density	flux_density	Ni2MnGa_sample_1	128	0.9855623999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_129_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_129_flux_density	flux_density	Ni2MnGa_sample_1	129	1.000047	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_130_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_130_flux_density	flux_density	Ni2MnGa_sample_1	130	1.0145315	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_131_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_131_flux_density	flux_density	Ni2MnGa_sample_1	131	1.029016	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_132_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_132_flux_density	flux_density	Ni2MnGa_sample_1	132	1.0435007	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_133_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_133_flux_density	flux_density	Ni2MnGa_sample_1	133	1.0579851999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_134_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_134_flux_density	flux_density	Ni2MnGa_sample_1	134	1.0724696999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_135_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_135_flux_density	flux_density	Ni2MnGa_sample_1	135	1.0869542000000003	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_136_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_136_flux_density	flux_density	Ni2MnGa_sample_1	136	1.1014389000000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_137_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_137_flux_density	flux_density	Ni2MnGa_sample_1	137	1.1159234	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_138_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_138_flux_density	flux_density	Ni2MnGa_sample_1	138	1.1304079000000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_139_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_139_flux_density	flux_density	Ni2MnGa_sample_1	139	1.1448924999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_140_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_140_flux_density	flux_density	Ni2MnGa_sample_1	140	1.1593771	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_141_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_141_flux_density	flux_density	Ni2MnGa_sample_1	141	1.1738616	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_142_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_142_flux_density	flux_density	Ni2MnGa_sample_1	142	1.1883461	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_143_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_143_flux_density	flux_density	Ni2MnGa_sample_1	143	1.2028308	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_144_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_144_flux_density	flux_density	Ni2MnGa_sample_1	144	1.2173153	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_145_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_145_flux_density	flux_density	Ni2MnGa_sample_1	145	1.2317998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_146_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_146_flux_density	flux_density	Ni2MnGa_sample_1	146	1.2462844	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_147_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_147_flux_density	flux_density	Ni2MnGa_sample_1	147	1.260769	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_148_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_148_flux_density	flux_density	Ni2MnGa_sample_1	148	1.2752535000000003	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_149_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_149_flux_density	flux_density	Ni2MnGa_sample_1	149	1.2897379999999996	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_150_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_150_flux_density	flux_density	Ni2MnGa_sample_1	150	1.3042226999999997	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_151_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_151_flux_density	flux_density	Ni2MnGa_sample_1	151	1.3187072	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_152_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_152_flux_density	flux_density	Ni2MnGa_sample_1	152	1.3331918	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_153_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_153_flux_density	flux_density	Ni2MnGa_sample_1	153	1.3476762999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_154_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_154_flux_density	flux_density	Ni2MnGa_sample_1	154	1.3621609	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_155_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_155_flux_density	flux_density	Ni2MnGa_sample_1	155	1.3766454	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_156_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_156_flux_density	flux_density	Ni2MnGa_sample_1	156	1.39113	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_157_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_157_flux_density	flux_density	Ni2MnGa_sample_1	157	1.4036959999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_158_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_158_flux_density	flux_density	Ni2MnGa_sample_1	158	1.4162620000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_159_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_159_flux_density	flux_density	Ni2MnGa_sample_1	159	1.428828	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_160_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_160_flux_density	flux_density	Ni2MnGa_sample_1	160	1.441394	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_161_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_161_flux_density	flux_density	Ni2MnGa_sample_1	161	1.4539600000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_162_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_162_flux_density	flux_density	Ni2MnGa_sample_1	162	1.466526	calc_flux_density	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_0magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1102	0	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_10000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1104	10000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_20000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1106	20000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_30000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1108	30000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_40000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1110	40000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_50000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1112	50000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_60000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1114	60000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_70000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1116	70000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_80000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1118	80000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_90000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1120	90000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_100000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1122	100000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_110000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1124	110000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_120000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1126	120000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_130000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1128	130000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_140000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1130	140000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_150000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1132	150000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_160000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1134	160000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_170000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1136	170000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_180000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1138	180000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_190000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1140	190000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_200000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1142	200000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_210000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1144	210000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_220000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1146	220000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_230000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1148	230000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_240000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1150	240000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_250000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1152	250000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_260000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1154	260000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_270000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1156	270000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_280000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1158	280000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_290000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1160	290000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_300000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1162	300000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_310000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1164	310000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_320000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1166	320000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_330000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1168	330000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_340000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1170	340000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_350000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1172	350000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_360000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1174	360000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_370000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1176	370000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_380000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1178	380000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_390000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1180	390000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_400000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1182	400000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_410000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1184	410000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_420000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1186	420000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_430000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1188	430000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_440000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1190	440000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_450000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1192	450000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_460000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1194	460000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_470000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1196	470000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_480000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1198	480000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_490000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1200	490000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_500000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1202	500000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_510000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1204	510000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_520000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1206	520000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_530000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1208	530000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_540000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1210	540000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_550000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1212	550000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_560000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1214	560000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_570000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1216	570000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_580000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1218	580000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_590000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1220	590000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_600000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1222	600000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_610000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1224	610000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_0	\N	magnetic_stress	Ni2MnGa_sample_1	1102	0	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_10000	\N	magnetic_stress	Ni2MnGa_sample_1	1104	23802.122937500026	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_20000	\N	magnetic_stress	Ni2MnGa_sample_1	1106	82526.78395000008	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_30000	\N	magnetic_stress	Ni2MnGa_sample_1	1108	151712.94220000017	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_40000	\N	magnetic_stress	Ni2MnGa_sample_1	1110	219581.26907500022	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_50000	\N	magnetic_stress	Ni2MnGa_sample_1	1112	286131.7667875003	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_60000	\N	magnetic_stress	Ni2MnGa_sample_1	1114	351364.4323875004	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_70000	\N	magnetic_stress	Ni2MnGa_sample_1	1116	415279.2666125005	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_80000	\N	magnetic_stress	Ni2MnGa_sample_1	1118	477876.2687250006	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_90000	\N	magnetic_stress	Ni2MnGa_sample_1	1120	539155.4387250006	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_100000	\N	magnetic_stress	Ni2MnGa_sample_1	1122	599116.7773500007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_110000	\N	magnetic_stress	Ni2MnGa_sample_1	1124	657760.2846000007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_120000	\N	magnetic_stress	Ni2MnGa_sample_1	1126	715085.9590000008	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_130000	\N	magnetic_stress	Ni2MnGa_sample_1	1128	771093.8005500009	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_140000	\N	magnetic_stress	Ni2MnGa_sample_1	1130	825783.814412501	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_150000	\N	magnetic_stress	Ni2MnGa_sample_1	1132	879155.9976375011	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_160000	\N	magnetic_stress	Ni2MnGa_sample_1	1134	931210.3472750011	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_170000	\N	magnetic_stress	Ni2MnGa_sample_1	1136	981946.864062501	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_180000	\N	magnetic_stress	Ni2MnGa_sample_1	1138	1031365.552425001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_190000	\N	magnetic_stress	Ni2MnGa_sample_1	1140	1079466.410150001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_200000	\N	magnetic_stress	Ni2MnGa_sample_1	1142	1126249.4350250012	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_210000	\N	magnetic_stress	Ni2MnGa_sample_1	1144	1171714.6285250012	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_220000	\N	magnetic_stress	Ni2MnGa_sample_1	1146	1215861.9877000013	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_230000	\N	magnetic_stress	Ni2MnGa_sample_1	1148	1258691.5184500013	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_240000	\N	magnetic_stress	Ni2MnGa_sample_1	1150	1300203.217825001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_250000	\N	magnetic_stress	Ni2MnGa_sample_1	1152	1340397.085825001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_260000	\N	magnetic_stress	Ni2MnGa_sample_1	1154	1379273.1202375009	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_270000	\N	magnetic_stress	Ni2MnGa_sample_1	1156	1416831.318850001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_280000	\N	magnetic_stress	Ni2MnGa_sample_1	1158	1453071.695675001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_290000	\N	magnetic_stress	Ni2MnGa_sample_1	1160	1487994.242600001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_300000	\N	magnetic_stress	Ni2MnGa_sample_1	1162	1521598.9463500008	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_310000	\N	magnetic_stress	Ni2MnGa_sample_1	1164	1553885.8216750007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_320000	\N	magnetic_stress	Ni2MnGa_sample_1	1166	1584854.8737375007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_330000	\N	magnetic_stress	Ni2MnGa_sample_1	1168	1614506.0870500007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_340000	\N	magnetic_stress	Ni2MnGa_sample_1	1170	1642839.457925001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_350000	\N	magnetic_stress	Ni2MnGa_sample_1	1172	1669855.004800001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_360000	\N	magnetic_stress	Ni2MnGa_sample_1	1174	1695552.724725001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_370000	\N	magnetic_stress	Ni2MnGa_sample_1	1176	1719932.6059000008	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_380000	\N	magnetic_stress	Ni2MnGa_sample_1	1178	1742994.658650001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_390000	\N	magnetic_stress	Ni2MnGa_sample_1	1180	1764738.887400001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_400000	\N	magnetic_stress	Ni2MnGa_sample_1	1182	1785165.2774000007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_410000	\N	magnetic_stress	Ni2MnGa_sample_1	1184	1804273.8212750005	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_420000	\N	magnetic_stress	Ni2MnGa_sample_1	1186	1822064.5411500004	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_430000	\N	magnetic_stress	Ni2MnGa_sample_1	1188	1838537.437025	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_440000	\N	magnetic_stress	Ni2MnGa_sample_1	1190	1853692.5015250004	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_450000	\N	magnetic_stress	Ni2MnGa_sample_1	1192	1867529.7346500005	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_460000	\N	magnetic_stress	Ni2MnGa_sample_1	1194	1880049.1364000007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_470000	\N	magnetic_stress	Ni2MnGa_sample_1	1196	1891250.7067750003	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_480000	\N	magnetic_stress	Ni2MnGa_sample_1	1198	1901134.4384	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_490000	\N	magnetic_stress	Ni2MnGa_sample_1	1200	1909700.3386499998	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_500000	\N	magnetic_stress	Ni2MnGa_sample_1	1202	1916948.4075250002	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_510000	\N	magnetic_stress	Ni2MnGa_sample_1	1204	1922878.645025	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_520000	\N	magnetic_stress	Ni2MnGa_sample_1	1206	1927491.05115	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_530000	\N	magnetic_stress	Ni2MnGa_sample_1	1208	1930785.6332750004	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_540000	\N	magnetic_stress	Ni2MnGa_sample_1	1210	1932762.3840250003	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_550000	\N	magnetic_stress	Ni2MnGa_sample_1	1212	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_560000	\N	magnetic_stress	Ni2MnGa_sample_1	1214	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_570000	\N	magnetic_stress	Ni2MnGa_sample_1	1216	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_580000	\N	magnetic_stress	Ni2MnGa_sample_1	1218	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_590000	\N	magnetic_stress	Ni2MnGa_sample_1	1220	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_600000	\N	magnetic_stress	Ni2MnGa_sample_1	1222	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_610000	\N	magnetic_stress	Ni2MnGa_sample_1	1224	1933421.2960250007	calc_magnetic_stress	4
cp_mod_calc_max_magn_stress_max_magnetic_stress_Ni2MnGa_sample_1_0	mod_calc_max_magn_stress_max_magnetic_stress_Ni2MnGa_sample_1	max_magnetic_stress	Ni2MnGa_sample_1	0	1933421.2960250007	calc_max_magn_stress	5
cp_mod_maximal_blocking_stress_load_calc_max_block_load_hold_stick_Ni2MnGa_sample_1_0_maximal_blocking_stress_hold	mod_calc_max_block_load_hold_maximal_blocking_stress_load_stick_Ni2MnGa_sample_1_0_maximal_blocking_stress_hold	maximal_blocking_stress_hold	stick_Ni2MnGa_sample_1	0	2233421.2960250005	calc_max_block_load_hold	6
cp_mod_maximal_blocking_stress_load_calc_max_block_load_hold_stick_Ni2MnGa_sample_1_0_maximal_blocking_stress_load	mod_calc_max_block_load_hold_maximal_blocking_stress_load_stick_Ni2MnGa_sample_1_0_maximal_blocking_stress_load	maximal_blocking_stress_load	stick_Ni2MnGa_sample_1	0	1633421.2960250007	calc_max_block_load_hold	6
\.


--
-- Data for Name: temp_view_results; Type: TABLE DATA; Schema: public; Owner: mena
--

COPY public.temp_view_results (concr_param_id, mod_meas_id, param_id, mat_sample_id, meas_time, value, mod_id, step) FROM stdin;
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_110_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_110_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	110	0	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_111_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_111_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	111	267550000	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_210_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_210_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	210	0	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_211_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_211_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	211	190775000	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_310_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_310_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	310	0	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v1_wire_actuator_NiTi#6_311_mechanic_stress	mod_1d_non_linear_mech_mod_v1_mechanic_stress_wire_actuator_NiTi#6_311_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	311	114000000	1d_non_linear_mech_mod_v1	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_112_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_112_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	112	471068000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_113_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_113_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	113	471818000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_114_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_114_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	114	472568000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_115_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_115_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	115	473318000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_116_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_116_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	116	474068000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_117_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_117_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	117	474818000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_118_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_118_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	118	475568000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_119_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_119_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	119	476318000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_120_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_120_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	120	477068000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_121_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_121_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	121	477818000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_212_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_212_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	212	335944000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_213_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_213_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	213	336694000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_214_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_214_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	214	337444000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_215_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_215_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	215	338194000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_216_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_216_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	216	338944000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_217_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_217_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	217	339694000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_218_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_218_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	218	340444000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_219_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_219_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	219	341194000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_220_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_220_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	220	341944000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_221_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_221_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	221	342694000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_312_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_312_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	312	200820000	1d_non_linear_mech_mod_v2	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_117_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_117_magnetization	magnetization	Ni2MnGa_sample_1	117	497514.1254177941	calc_magnetization	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_313_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_313_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	313	201570000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_314_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_314_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	314	202320000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_315_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_315_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	315	203070000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_316_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_316_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	316	203820000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_317_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_317_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	317	204570000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_318_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_318_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	318	205320000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_319_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_319_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	319	206070000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_320_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_320_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	320	206820000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v2_wire_actuator_NiTi#6_321_mechanic_stress	mod_1d_non_linear_mech_mod_v2_mechanic_stress_wire_actuator_NiTi#6_321_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	321	207570000	1d_non_linear_mech_mod_v2	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v3_wire_actuator_NiTi#6_122_mechanic_stress	mod_1d_non_linear_mech_mod_v3_mechanic_stress_wire_actuator_NiTi#6_122_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	122	525745000	1d_non_linear_mech_mod_v3	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v3_wire_actuator_NiTi#6_222_mechanic_stress	mod_1d_non_linear_mech_mod_v3_mechanic_stress_wire_actuator_NiTi#6_222_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	222	390621000	1d_non_linear_mech_mod_v3	1
cp_mod_mechanic_stress_1d_non_linear_mech_mod_v3_wire_actuator_NiTi#6_322_mechanic_stress	mod_1d_non_linear_mech_mod_v3_mechanic_stress_wire_actuator_NiTi#6_322_mechanic_stress	mechanic_stress	wire_actuator_NiTi#6	322	255496999.99999997	1d_non_linear_mech_mod_v3	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_1_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_1_magnetization	magnetization	Ni2MnGa_sample_1	1	0	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_10_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_10_magnetization	magnetization	Ni2MnGa_sample_1	10	159765.68518223776	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_101_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_101_magnetization	magnetization	Ni2MnGa_sample_1	101	0	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_102_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_102_magnetization	magnetization	Ni2MnGa_sample_1	102	336629.15804551967	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_103_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_103_magnetization	magnetization	Ni2MnGa_sample_1	103	476139.12939678505	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_104_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_104_magnetization	magnetization	Ni2MnGa_sample_1	104	477665.89208976604	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_105_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_105_magnetization	magnetization	Ni2MnGa_sample_1	105	479192.69457265636	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_106_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_106_magnetization	magnetization	Ni2MnGa_sample_1	106	480719.48113958305	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_107_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_107_magnetization	magnetization	Ni2MnGa_sample_1	107	482246.2597485278	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_108_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_108_magnetization	magnetization	Ni2MnGa_sample_1	108	483773.0622314181	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_109_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_109_magnetization	magnetization	Ni2MnGa_sample_1	109	485299.82492439914	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_11_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_11_magnetization	magnetization	Ni2MnGa_sample_1	11	168402.47493235715	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_110_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_110_magnetization	magnetization	Ni2MnGa_sample_1	110	486826.62740728946	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_111_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_111_magnetization	magnetization	Ni2MnGa_sample_1	111	488353.3901002706	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_112_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_112_magnetization	magnetization	Ni2MnGa_sample_1	112	489880.19258316094	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_113_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_113_magnetization	magnetization	Ni2MnGa_sample_1	113	491406.9711921057	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_114_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_114_magnetization	magnetization	Ni2MnGa_sample_1	114	492933.7577590323	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_115_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_115_magnetization	magnetization	Ni2MnGa_sample_1	115	494460.5681999045	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_116_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_116_magnetization	magnetization	Ni2MnGa_sample_1	116	495987.3229349038	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_118_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_118_magnetization	magnetization	Ni2MnGa_sample_1	118	499040.90402673883	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_119_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_119_magnetization	magnetization	Ni2MnGa_sample_1	119	500567.69059366547	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_12_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_12_magnetization	magnetization	Ni2MnGa_sample_1	12	177039.2885564221	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_120_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_120_magnetization	magnetization	Ni2MnGa_sample_1	120	502094.4612446284	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_121_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_121_magnetization	magnetization	Ni2MnGa_sample_1	121	503621.2796434824	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_122_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_122_magnetization	magnetization	Ni2MnGa_sample_1	122	505148.06621040904	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_123_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_123_magnetization	magnetization	Ni2MnGa_sample_1	123	506674.836861372	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_124_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_124_magnetization	magnetization	Ni2MnGa_sample_1	124	508201.6313862804	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_125_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_125_magnetization	magnetization	Ni2MnGa_sample_1	125	509728.3940792616	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_126_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_126_magnetization	magnetization	Ni2MnGa_sample_1	126	511255.2124781156	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_127_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_127_magnetization	magnetization	Ni2MnGa_sample_1	127	512781.95129715104	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_128_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_128_magnetization	magnetization	Ni2MnGa_sample_1	128	514308.76969600504	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_129_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_129_magnetization	magnetization	Ni2MnGa_sample_1	129	515835.58809485915	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_13_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_13_magnetization	magnetization	Ni2MnGa_sample_1	13	185676.1101384689	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_130_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_130_magnetization	magnetization	Ni2MnGa_sample_1	130	517362.3269138946	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_131_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_131_magnetization	magnetization	Ni2MnGa_sample_1	131	518889.0657329301	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_132_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_132_magnetization	magnetization	Ni2MnGa_sample_1	132	520415.96371160285	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_133_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_133_magnetization	magnetization	Ni2MnGa_sample_1	133	521942.7025306383	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_134_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_134_magnetization	magnetization	Ni2MnGa_sample_1	134	523469.4413496738	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_135_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_135_magnetization	magnetization	Ni2MnGa_sample_1	135	524996.1801687094	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_136_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_136_magnetization	magnetization	Ni2MnGa_sample_1	136	526523.0781473819	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_137_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_137_magnetization	magnetization	Ni2MnGa_sample_1	137	528049.8169664174	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_138_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_138_magnetization	magnetization	Ni2MnGa_sample_1	138	529576.5557854528	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_139_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_139_magnetization	magnetization	Ni2MnGa_sample_1	139	531103.374184307	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_14_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_14_magnetization	magnetization	Ni2MnGa_sample_1	14	194312.90784657013	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_140_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_140_magnetization	magnetization	Ni2MnGa_sample_1	140	532630.1925831609	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_141_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_141_magnetization	magnetization	Ni2MnGa_sample_1	141	534156.9314021964	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_142_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_142_magnetization	magnetization	Ni2MnGa_sample_1	142	535683.6702212319	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_143_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_143_magnetization	magnetization	Ni2MnGa_sample_1	143	537210.5681999046	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_144_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_144_magnetization	magnetization	Ni2MnGa_sample_1	144	538737.3070189401	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_145_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_145_magnetization	magnetization	Ni2MnGa_sample_1	145	540264.0458379755	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_146_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_146_magnetization	magnetization	Ni2MnGa_sample_1	146	541790.8642368296	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_147_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_147_magnetization	magnetization	Ni2MnGa_sample_1	147	543317.6826356837	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_148_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_148_magnetization	magnetization	Ni2MnGa_sample_1	148	544844.4214547192	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_149_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_149_magnetization	magnetization	Ni2MnGa_sample_1	149	546371.1602737546	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_15_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_15_magnetization	magnetization	Ni2MnGa_sample_1	15	202949.70555467135	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_150_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_150_magnetization	magnetization	Ni2MnGa_sample_1	150	547898.0582524271	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_151_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_151_magnetization	magnetization	Ni2MnGa_sample_1	151	549424.7970714627	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_152_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_152_magnetization	magnetization	Ni2MnGa_sample_1	152	550951.6154703167	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_153_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_153_magnetization	magnetization	Ni2MnGa_sample_1	153	552478.3542893522	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_154_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_154_magnetization	magnetization	Ni2MnGa_sample_1	154	554005.1726882062	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_155_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_155_magnetization	magnetization	Ni2MnGa_sample_1	155	555531.9115072417	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_156_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_156_magnetization	magnetization	Ni2MnGa_sample_1	156	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_157_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_157_magnetization	magnetization	Ni2MnGa_sample_1	157	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_158_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_158_magnetization	magnetization	Ni2MnGa_sample_1	158	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_159_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_159_magnetization	magnetization	Ni2MnGa_sample_1	159	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_16_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_16_magnetization	magnetization	Ni2MnGa_sample_1	16	211586.50326277257	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_160_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_160_magnetization	magnetization	Ni2MnGa_sample_1	160	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_161_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_161_magnetization	magnetization	Ni2MnGa_sample_1	161	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_162_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_162_magnetization	magnetization	Ni2MnGa_sample_1	162	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_17_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_17_magnetization	magnetization	Ni2MnGa_sample_1	17	220223.32484481938	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_18_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_18_magnetization	magnetization	Ni2MnGa_sample_1	18	228860.13846888428	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_19_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_19_magnetization	magnetization	Ni2MnGa_sample_1	19	237496.8963870763	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_2_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_2_magnetization	magnetization	Ni2MnGa_sample_1	2	79792.73436256566	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_20_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_20_magnetization	magnetization	Ni2MnGa_sample_1	20	246133.72592710488	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_21_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_21_magnetization	magnetization	Ni2MnGa_sample_1	21	254770.5395511698	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_22_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_22_magnetization	magnetization	Ni2MnGa_sample_1	22	263407.3690911985	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_23_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_23_magnetization	magnetization	Ni2MnGa_sample_1	23	272044.1667992997	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_24_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_24_magnetization	magnetization	Ni2MnGa_sample_1	24	280680.9406334554	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_25_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_25_magnetization	magnetization	Ni2MnGa_sample_1	25	289317.7622155022	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_26_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_26_magnetization	magnetization	Ni2MnGa_sample_1	26	297954.55992360343	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_27_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_27_magnetization	magnetization	Ni2MnGa_sample_1	27	306591.38150565024	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_28_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_28_magnetization	magnetization	Ni2MnGa_sample_1	28	315228.20308769704	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_29_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_29_magnetization	magnetization	Ni2MnGa_sample_1	29	323864.9530479071	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_3_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_3_magnetization	magnetization	Ni2MnGa_sample_1	3	99308.04551965621	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_30_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_30_magnetization	magnetization	Ni2MnGa_sample_1	30	332501.7825879357	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_31_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_31_magnetization	magnetization	Ni2MnGa_sample_1	31	341138.5962120007	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_32_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_32_magnetization	magnetization	Ni2MnGa_sample_1	32	349775.4257520293	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_33_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_33_magnetization	magnetization	Ni2MnGa_sample_1	33	358412.18367022125	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_34_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_34_magnetization	magnetization	Ni2MnGa_sample_1	34	367049.02116823173	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_35_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_35_magnetization	magnetization	Ni2MnGa_sample_1	35	375685.81887633295	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_36_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_36_magnetization	magnetization	Ni2MnGa_sample_1	36	384322.6165844342	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_37_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_37_magnetization	magnetization	Ni2MnGa_sample_1	37	392959.4461244629	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_38_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_38_magnetization	magnetization	Ni2MnGa_sample_1	38	401596.2120006366	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_39_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_39_magnetization	magnetization	Ni2MnGa_sample_1	39	410233.00970873784	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_4_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_4_magnetization	magnetization	Ni2MnGa_sample_1	4	107944.8511857393	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_40_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_40_magnetization	magnetization	Ni2MnGa_sample_1	40	418869.80741683906	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_41_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_41_magnetization	magnetization	Ni2MnGa_sample_1	41	427506.6847047589	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_42_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_42_magnetization	magnetization	Ni2MnGa_sample_1	42	436143.48241286003	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_43_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_43_magnetization	magnetization	Ni2MnGa_sample_1	43	444780.28012096137	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_44_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_44_magnetization	magnetization	Ni2MnGa_sample_1	44	453417.0778290626	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_45_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_45_magnetization	magnetization	Ni2MnGa_sample_1	45	462053.7959573453	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_46_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_46_magnetization	magnetization	Ni2MnGa_sample_1	46	470690.67324526503	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_47_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_47_magnetization	magnetization	Ni2MnGa_sample_1	47	479327.47095336637	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_48_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_48_magnetization	magnetization	Ni2MnGa_sample_1	48	487964.2686614675	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_49_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_49_magnetization	magnetization	Ni2MnGa_sample_1	49	496601.0663695687	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_5_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_5_magnetization	magnetization	Ni2MnGa_sample_1	5	116581.64889384052	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_50_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_50_magnetization	magnetization	Ni2MnGa_sample_1	50	505237.94365748845	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_51_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_51_magnetization	magnetization	Ni2MnGa_sample_1	51	513874.7413655898	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_52_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_52_magnetization	magnetization	Ni2MnGa_sample_1	52	522511.539073691	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_53_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_53_magnetization	magnetization	Ni2MnGa_sample_1	53	531148.3367817921	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_54_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_54_magnetization	magnetization	Ni2MnGa_sample_1	54	539785.0549100749	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_55_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_55_magnetization	magnetization	Ni2MnGa_sample_1	55	548421.9321979946	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_56_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_56_magnetization	magnetization	Ni2MnGa_sample_1	56	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_57_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_57_magnetization	magnetization	Ni2MnGa_sample_1	57	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_58_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_58_magnetization	magnetization	Ni2MnGa_sample_1	58	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_59_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_59_magnetization	magnetization	Ni2MnGa_sample_1	59	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_6_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_6_magnetization	magnetization	Ni2MnGa_sample_1	6	125218.45455992362	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_60_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_60_magnetization	magnetization	Ni2MnGa_sample_1	60	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_61_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_61_magnetization	magnetization	Ni2MnGa_sample_1	61	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_62_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_62_magnetization	magnetization	Ni2MnGa_sample_1	62	557058.7299060959	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_7_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_7_magnetization	magnetization	Ni2MnGa_sample_1	7	133855.2602260067	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_8_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_8_magnetization	magnetization	Ni2MnGa_sample_1	8	142492.07385007164	calc_magnetization	1
cp_mod_magnetization_calc_magnetization_Ni2MnGa_sample_1_9_magnetization	mod_calc_magnetization_magnetization_Ni2MnGa_sample_1_9_magnetization	magnetization	Ni2MnGa_sample_1	9	151128.87155817286	calc_magnetization	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_131_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_131_martensite_content	martensite_content	wire_actuator_1	131	1	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_132_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_132_martensite_content	martensite_content	wire_actuator_1	132	0.9427280128266052	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_133_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_133_martensite_content	martensite_content	wire_actuator_1	133	0.43973165987233864	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_134_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_134_martensite_content	martensite_content	wire_actuator_1	134	0.16843867087960235	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_135_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_135_martensite_content	martensite_content	wire_actuator_1	135	0.014529091286973939	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_136_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_136_martensite_content	martensite_content	wire_actuator_1	136	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_137_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_137_martensite_content	martensite_content	wire_actuator_1	137	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_138_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_138_martensite_content	martensite_content	wire_actuator_1	138	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_139_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_139_martensite_content	martensite_content	wire_actuator_1	139	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_140_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_140_martensite_content	martensite_content	wire_actuator_1	140	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_141_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_141_martensite_content	martensite_content	wire_actuator_1	141	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_142_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_142_martensite_content	martensite_content	wire_actuator_1	142	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_asc_wire_actuator_1_143_martensite_content	mod_calc_martensite_asc_martensite_content_wire_actuator_1_143_martensite_content	martensite_content	wire_actuator_1	143	0	calc_martensite_asc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_111_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_111_martensite_content	martensite_content	wire_actuator_1	111	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_112_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_112_martensite_content	martensite_content	wire_actuator_1	112	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_113_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_113_martensite_content	martensite_content	wire_actuator_1	113	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_114_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_114_martensite_content	martensite_content	wire_actuator_1	114	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_115_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_115_martensite_content	martensite_content	wire_actuator_1	115	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_116_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_116_martensite_content	martensite_content	wire_actuator_1	116	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_117_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_117_martensite_content	martensite_content	wire_actuator_1	117	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_118_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_118_martensite_content	martensite_content	wire_actuator_1	118	1	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_119_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_119_martensite_content	martensite_content	wire_actuator_1	119	0.9567727288213013	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_120_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_120_martensite_content	martensite_content	wire_actuator_1	120	0.7191855733945405	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_121_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_121_martensite_content	martensite_content	wire_actuator_1	121	0.37903905220016804	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_122_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_122_martensite_content	martensite_content	wire_actuator_1	122	0.019369152030841108	calc_martensite_desc	1
cp_mod_martensite_content_calc_martensite_desc_wire_actuator_1_123_martensite_content	mod_calc_martensite_desc_martensite_content_wire_actuator_1_123_martensite_content	martensite_content	wire_actuator_1	123	0	calc_martensite_desc	1
cp_mod_maximum_strain_mag_calc_maximum_strain_mag_Ni2MnGa_sample_1_0_maximum_strain_mag	mod_calc_maximum_strain_mag_maximum_strain_mag_Ni2MnGa_sample_1_0_maximum_strain_mag	maximum_strain_mag	Ni2MnGa_sample_1	0	0.06779661016949146	calc_maximum_strain_mag	1
cp_mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_1	mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_1	Yeoh_vector_1	sample1_elastomer1	0	194409.9644002127	calc_yeoh_minR	1
cp_mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_2	mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_2	Yeoh_vector_2	sample1_elastomer1	0	-2275.4747645269717	calc_yeoh_minR	1
cp_mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_3	mod_calc_yeoh_minR_sample1_elastomer1Yeoh_vector_3	Yeoh_vector_3	sample1_elastomer1	0	527.8018381989875	calc_yeoh_minR	1
cp_mod_calc_young_modulus_minR_sample1_elastomer1young_modulus	mod_calc_young_modulus_minR_sample1_elastomer1young_modulus	young_modulus	sample1_elastomer1	0	652412.0602756739	calc_young_modulus_minR	1
cp_mod_calc_young_modulus_Neo_Hookean_minR_sample1_elastomer1young_modulus_neo_Hookean	mod_calc_young_modulus_Neo_Hookean_minR_sample1_elastomer1young_modulus_neo_Hookean	young_modulus_neo_Hookean	sample1_elastomer1	0	1178216.3521274924	calc_young_modulus_Neo_Hookean_minR	1
cp_mod_maximal_blocking_stress_electrostatic_pressure_model_elastomer_1_0_maximal_blocking_stress	mod_electrostatic_pressure_model_maximal_blocking_stress_elastomer_1_0_maximal_blocking_stress	maximal_blocking_stress	elastomer_1	0	277689.7728	electrostatic_pressure_model	1
cp_mod_linear_interpol_initial_stress_initial_stress_linear_interpolation_initial_stress_sample1_elastomer1_0	mod_linear_interpolation_initial_stress_initial_stress_sample1_elastomer1_0	initial_stress	sample1_elastomer1	0	59400	linear_interpolation_initial_stress	1
cp_mod_A10_matrix_model_1_M420_0_1	mod_matrix_model_1_A10_M420_0	A10	M420	0	120197869377.449	matrix_model_1	1
cp_mod_A10_matrix_model_1_M420_0_2	mod_matrix_model_1_A10_M420_0	A10	M420	0	72804504448.5392	matrix_model_1	1
cp_mod_A10_matrix_model_1_M420_0_3	mod_matrix_model_1_A10_M420_0	A10	M420	0	67086386624.0066	matrix_model_1	1
cp_mod_A10_matrix_model_1_M420_0_4	mod_matrix_model_1_A10_M420_0	A10	M420	0	100113530808.133	matrix_model_1	1
cp_mod_A10_matrix_model_1_M420_0_5	mod_matrix_model_1_A10_M420_0	A10	M420	0	22222222222.2222	matrix_model_1	1
cp_mod_B10_matrix_model_2_M420_0_1	mod_matrix_model_2_B10_M420_0	B10	M420	0	-7.06471256063577	matrix_model_2	1
cp_mod_B10_matrix_model_2_M420_0_2	mod_matrix_model_2_B10_M420_0	B10	M420	0	14.0726597172051	matrix_model_2	1
cp_mod_B10_matrix_model_2_M420_0_3	mod_matrix_model_2_B10_M420_0	B10	M420	0	11.6666666666667	matrix_model_2	1
cp_mod_C10_matrix_model_3_M420_0_1	mod_matrix_model_3_C10_M420_0	C10	M420	0	8.0414e-09	matrix_model_3	1
cp_mod_C10_matrix_model_3_M420_0_2	mod_matrix_model_3_C10_M420_0	C10	M420	0	6.90989778098875e-09	matrix_model_3	1
cp_mod_A01_matrix_model_4_M420_0_1	mod_matrix_model_4_A01_M420_0	A01	M420	0	1.35929071606054e-11	matrix_model_4	1
cp_mod_A01_matrix_model_4_M420_0_2	mod_matrix_model_4_A01_M420_0	A01	M420	0	-7.50709283939462e-12	matrix_model_4	1
cp_mod_A01_matrix_model_4_M420_0_3	mod_matrix_model_4_A01_M420_0	A01	M420	0	-2.49051276259318e-12	matrix_model_4	1
cp_mod_A01_matrix_model_4_M420_0_4	mod_matrix_model_4_A01_M420_0	A01	M420	0	9.80395019200361e-12	matrix_model_4	1
cp_mod_A01_matrix_model_4_M420_0_5	mod_matrix_model_4_A01_M420_0	A01	M420	0	2.55437514117913e-11	matrix_model_4	1
cp_mod_B01_matrix_model_5_M420_0_1	mod_matrix_model_5_B01_M420_0	B01	M420	0	-0.0112943302462164	matrix_model_5	1
cp_mod_B01_matrix_model_5_M420_0_2	mod_matrix_model_5_B01_M420_0	B01	M420	0	0.0250592952337926	matrix_model_5	1
cp_mod_B01_matrix_model_5_M420_0_3	mod_matrix_model_5_B01_M420_0	B01	M420	0	0.0370595211203976	matrix_model_5	1
cp_mod_C01_matrix_model_6_M420_0_1	mod_matrix_model_6_C01_M420_0	C01	M420	0	70589564.0388525	matrix_model_6	1
cp_mod_C01_matrix_model_6_M420_0_2	mod_matrix_model_6_C01_M420_0	C01	M420	0	70589564.0388525	matrix_model_6	1
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_1_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_1_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	1	0	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_10_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_10_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	10	0.20076155999999998	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_101_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_101_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	101	0	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_102_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_102_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	102	0.4230082	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_103_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_103_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	103	0.5983164300000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_104_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_104_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	104	0.60023496	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_105_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_105_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	105	0.60215354	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_106_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_106_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	106	0.6040721	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_107_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_107_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	107	0.60599065	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_108_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_108_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	108	0.60790923	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_109_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_109_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	109	0.6098277599999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_11_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_11_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	11	0.21161454999999998	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_110_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_110_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	110	0.6117463399999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_111_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_111_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	111	0.61366487	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_112_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_112_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	112	0.61558345	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_113_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_113_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	113	0.617502	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_114_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_114_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	114	0.61942056	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_115_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_115_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	115	0.62133915	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_116_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_116_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	116	0.6232576700000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_117_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_117_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	117	0.6251762500000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_118_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_118_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	118	0.6270948	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_119_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_119_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	119	0.62901336	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_12_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_12_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	12	0.22246757	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_120_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_120_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	120	0.6309319000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_121_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_121_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	121	0.6328505	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_122_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_122_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	122	0.63476906	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_123_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_123_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	123	0.6366876	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_124_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_124_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	124	0.63860617	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_125_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_125_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	125	0.6405247000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_126_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_126_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	126	0.6424433	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_127_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_127_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	127	0.6443618	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_128_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_128_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	128	0.6462803999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_129_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_129_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	129	0.648199	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_13_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_13_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	13	0.23332060000000002	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_130_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_130_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	130	0.6501175	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_131_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_131_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	131	0.652036	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_132_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_132_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	132	0.6539547000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_133_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_133_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	133	0.6558732	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_134_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_134_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	134	0.6577917000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_135_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_135_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	135	0.6597102000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_136_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_136_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	136	0.6616289000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_137_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_137_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	137	0.6635474	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_138_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_138_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	138	0.6654659	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_139_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_139_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	139	0.6673845	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_14_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_14_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	14	0.24417360000000002	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_140_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_140_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	140	0.6693031	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_141_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_141_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	141	0.6712216	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_142_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_142_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	142	0.6731400999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_143_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_143_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	143	0.6750588000000002	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_144_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_144_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	144	0.6769773000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_145_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_145_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	145	0.6788957999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_146_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_146_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	146	0.6808144	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_147_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_147_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	147	0.6827330000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_148_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_148_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	148	0.6846515000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_149_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_149_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	149	0.68657	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_15_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_15_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	15	0.2550266	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_150_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_150_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	150	0.6884887	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_151_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_151_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	151	0.6904072	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_152_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_152_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	152	0.6923258	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_153_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_153_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	153	0.6942442999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_154_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_154_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	154	0.6961628999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_155_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_155_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	155	0.6980813999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_156_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_156_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	156	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_157_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_157_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	157	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_158_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_158_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	158	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_159_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_159_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	159	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_16_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_16_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	16	0.2658796	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_160_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_160_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	160	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_161_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_161_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	161	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_162_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_162_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	162	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_17_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_17_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	17	0.27673263000000003	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_18_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_18_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	18	0.28758564999999997	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_19_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_19_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	19	0.29843860000000005	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_2_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_2_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	2	0.10026755000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_20_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_20_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	20	0.30929164	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_21_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_21_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	21	0.32014465999999997	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_22_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_22_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	22	0.3309977	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_23_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_23_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	23	0.3418507	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_24_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_24_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	24	0.3527036700000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_25_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_25_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	25	0.36355670000000007	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_26_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_26_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	26	0.37440970000000007	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_27_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_27_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	27	0.38526273000000005	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_28_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_28_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	28	0.3961157600000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_29_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_29_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	29	0.40696870000000007	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_3_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_3_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	3	0.12479048999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_30_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_30_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	30	0.41782173999999994	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_31_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_31_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	31	0.42867476000000004	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_32_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_32_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	32	0.43952779999999997	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_33_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_33_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	33	0.45038075	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_34_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_34_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	34	0.46123379999999997	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_35_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_35_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	35	0.4720868	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_36_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_36_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	36	0.4829398	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_37_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_37_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	37	0.49379284	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_38_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_38_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	38	0.5046457999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_39_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_39_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	39	0.5154987999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_4_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_4_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	4	0.1356435	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_40_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_40_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	40	0.5263517999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_41_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_41_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	41	0.5372049000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_42_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_42_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	42	0.5480578999999999	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_43_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_43_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	43	0.5589109	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_44_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_44_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	44	0.5697639	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_45_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_45_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	45	0.5806168	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_46_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_46_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	46	0.5914699	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_47_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_47_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	47	0.6023229000000002	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_48_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_48_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	48	0.6131759	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_49_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_49_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	49	0.6240289	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_5_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_5_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	5	0.1464965	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_50_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_50_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	50	0.634882	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_51_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_51_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	51	0.6457350000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_52_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_52_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	52	0.6565880000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_53_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_53_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	53	0.667441	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_54_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_54_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	54	0.6782939000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_55_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_55_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	55	0.689147	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_56_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_56_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	56	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_57_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_57_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	57	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_58_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_58_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	58	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_59_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_59_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	59	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_6_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_6_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	6	0.15734951	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_60_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_60_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	60	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_61_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_61_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	61	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_62_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_62_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	62	0.7000000000000001	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_7_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_7_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	7	0.16820252000000002	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_8_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_8_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	8	0.17905554	calc_mag_polarization	2
cp_mod_magnetic_polarization_calc_mag_polarization_Ni2MnGa_sample_1_9_magnetic_polarization	mod_calc_mag_polarization_magnetic_polarization_Ni2MnGa_sample_1_9_magnetic_polarization	magnetic_polarization	Ni2MnGa_sample_1	9	0.18990854000000001	calc_mag_polarization	2
cp_mod_young_modulus_neo_Hookean_calc_neo_hookean_initial_slope_sample1_elastomer1_0_young_modulus_neo_Hookean	mod_calc_neo_hookean_initial_slope_young_modulus_neo_Hookean_sample1_elastomer1_0_young_modulus_neo_Hookean	young_modulus_neo_Hookean	sample1_elastomer1	0	1187057.8905630445	calc_neo_hookean_initial_slope	2
cp_mod_young_modulus_calc_young_modulus_initial_slope_sample1_elastomer1_0_young_modulus	mod_calc_young_modulus_initial_slope_young_modulus_sample1_elastomer1_0_young_modulus	young_modulus	sample1_elastomer1	0	1188000	calc_young_modulus_initial_slope	2
cp_mod_A00_matrix_model_1_M420_0_1	mod_matrix_model_1_A00_M420_0	A00	M420	0	1.54e-11	matrix_model_1	2
cp_mod_A00_matrix_model_1_M420_0_2	mod_matrix_model_1_A00_M420_0	A00	M420	0	-5.70000000000003e-12	matrix_model_1	2
cp_mod_A00_matrix_model_1_M420_0_3	mod_matrix_model_1_A00_M420_0	A00	M420	0	-6.49999999999999e-12	matrix_model_1	2
cp_mod_A00_matrix_model_1_M420_0_4	mod_matrix_model_1_A00_M420_0	A00	M420	0	1.87e-11	matrix_model_1	2
cp_mod_A00_matrix_model_1_M420_0_5	mod_matrix_model_1_A00_M420_0	A00	M420	0	4.5e-11	matrix_model_1	2
cp_mod_A11_matrix_model_1_M420_0_1	mod_matrix_model_1_A11_M420_0	A11	M420	0	127420865309.707	matrix_model_1	2
cp_mod_A11_matrix_model_1_M420_0_2	mod_matrix_model_1_A11_M420_0	A11	M420	0	80027500380.7973	matrix_model_1	2
cp_mod_A11_matrix_model_1_M420_0_3	mod_matrix_model_1_A11_M420_0	A11	M420	0	52698431980.2742	matrix_model_1	2
cp_mod_A11_matrix_model_1_M420_0_4	mod_matrix_model_1_A11_M420_0	A11	M420	0	128773831986.699	matrix_model_1	2
cp_mod_A11_matrix_model_1_M420_0_5	mod_matrix_model_1_A11_M420_0	A11	M420	0	39148517532.8784	matrix_model_1	2
cp_mod_B00_matrix_model_2_M420_0_1	mod_matrix_model_2_B00_M420_0	B00	M420	0	-1.6e-10	matrix_model_2	2
cp_mod_B00_matrix_model_2_M420_0_2	mod_matrix_model_2_B00_M420_0	B00	M420	0	3.55e-10	matrix_model_2	2
cp_mod_B00_matrix_model_2_M420_0_3	mod_matrix_model_2_B00_M420_0	B00	M420	0	5.25000000000002e-10	matrix_model_2	2
cp_mod_B11_matrix_model_2_M420_0_1	mod_matrix_model_2_B11_M420_0	B11	M420	0	-964471223.088693	matrix_model_2	2
cp_mod_B11_matrix_model_2_M420_0_2	mod_matrix_model_2_B11_M420_0	B11	M420	0	2059980296.10146	matrix_model_2	2
cp_mod_B11_matrix_model_2_M420_0_3	mod_matrix_model_2_B11_M420_0	B11	M420	0	1450825312.34196	matrix_model_2	2
cp_mod_C00_matrix_model_3_M420_0_1	mod_matrix_model_3_C00_M420_0	C00	M420	0	1.41664e-08	matrix_model_3	2
cp_mod_C00_matrix_model_3_M420_0_2	mod_matrix_model_3_C00_M420_0	C00	M420	0	1.41664e-08	matrix_model_3	2
cp_mod_C11_matrix_model_3_M420_0_1	mod_matrix_model_3_C11_M420_0	C11	M420	0	124356455.343597	matrix_model_3	2
cp_mod_C11_matrix_model_3_M420_0_2	mod_matrix_model_3_C11_M420_0	C11	M420	0	144719941.118566	matrix_model_3	2
cp_mod_A00_matrix_model_4_M420_0_1	mod_matrix_model_4_A00_M420_0	A00	M420	0	1.53070415606054e-11	matrix_model_4	2
cp_mod_A00_matrix_model_4_M420_0_2	mod_matrix_model_4_A00_M420_0	A00	M420	0	-5.79295843939462e-12	matrix_model_4	2
cp_mod_A00_matrix_model_4_M420_0_3	mod_matrix_model_4_A00_M420_0	A00	M420	0	-6.38627276259318e-12	matrix_model_4	2
cp_mod_A00_matrix_model_4_M420_0_4	mod_matrix_model_4_A00_M420_0	A00	M420	0	1.86579501920036e-11	matrix_model_4	2
cp_mod_A00_matrix_model_4_M420_0_5	mod_matrix_model_4_A00_M420_0	A00	M420	0	4.50000000000001e-11	matrix_model_4	2
cp_mod_A11_matrix_model_4_M420_0_1	mod_matrix_model_4_A11_M420_0	A11	M420	0	127420865309.708	matrix_model_4	2
cp_mod_A11_matrix_model_4_M420_0_2	mod_matrix_model_4_A11_M420_0	A11	M420	0	80027500380.7979	matrix_model_4	2
cp_mod_A11_matrix_model_4_M420_0_3	mod_matrix_model_4_A11_M420_0	A11	M420	0	52698431980.2744	matrix_model_4	2
cp_mod_A11_matrix_model_4_M420_0_4	mod_matrix_model_4_A11_M420_0	A11	M420	0	128773831986.7	matrix_model_4	2
cp_mod_A11_matrix_model_4_M420_0_5	mod_matrix_model_4_A11_M420_0	A11	M420	0	39148517532.8785	matrix_model_4	2
cp_mod_B00_matrix_model_5_M420_0_1	mod_matrix_model_5_B00_M420_0	B00	M420	0	-1.6e-10	matrix_model_5	2
cp_mod_B00_matrix_model_5_M420_0_2	mod_matrix_model_5_B00_M420_0	B00	M420	0	3.5416e-10	matrix_model_5	2
cp_mod_B00_matrix_model_5_M420_0_3	mod_matrix_model_5_B00_M420_0	B00	M420	0	5.25000000000001e-10	matrix_model_5	2
cp_mod_B11_matrix_model_5_M420_0_1	mod_matrix_model_5_B11_M420_0	B11	M420	0	-1022404785.7948	matrix_model_5	2
cp_mod_B11_matrix_model_5_M420_0_2	mod_matrix_model_5_B11_M420_0	B11	M420	0	2036594485.65553	matrix_model_5	2
cp_mod_B11_matrix_model_5_M420_0_3	mod_matrix_model_5_B11_M420_0	B11	M420	0	1450825312.34197	matrix_model_5	2
cp_mod_C00_matrix_model_6_M420_0_1	mod_matrix_model_6_C00_M420_0	C00	M420	0	1.41664e-08	matrix_model_6	2
cp_mod_C00_matrix_model_6_M420_0_2	mod_matrix_model_6_C00_M420_0	C00	M420	0	1.41664e-08	matrix_model_6	2
cp_mod_C11_matrix_model_6_M420_0_1	mod_matrix_model_6_C11_M420_0	C11	M420	0	124356455.343597	matrix_model_6	2
cp_mod_C11_matrix_model_6_M420_0_2	mod_matrix_model_6_C11_M420_0	C11	M420	0	144719941.118566	matrix_model_6	2
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_6	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	6	5858.72795	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_7	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	7	7486.4881000000005	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_1	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	1	0	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_2	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	2	501.33775	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_3	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	3	1626.62795	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_4	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	4	2928.7979	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_5	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	5	4339.4979	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_8	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	8	9222.778400000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_9	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	9	11067.598800000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_10	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	10	13020.949300000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_11	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	11	15082.829850000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_12	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	12	17253.24045	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_13	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	13	19532.1813	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_14	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	14	21919.6523	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_15	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	15	24415.6533	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_16	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	16	27020.1843	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_17	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	17	29733.245450000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_18	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	18	32554.836850000003	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_19	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	19	35484.9581	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_20	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	20	38523.609300000004	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_21	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	21	41670.7908	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_22	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	22	44926.5026	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_23	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	23	48290.7446	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_24	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	24	51763.516449999996	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_25	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	25	55344.8183	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_26	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	26	59034.6503	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_27	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	27	62833.01245	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_28	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	28	66739.90490000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_29	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	29	70755.32720000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_30	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	30	74879.27940000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_31	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	31	79111.76190000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_32	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	32	83452.77470000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_33	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	33	87902.31745	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_34	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	34	92460.3902	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_35	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	35	97126.9932	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_36	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	36	101902.1262	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_37	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	37	106785.7894	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_38	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	38	111777.98259999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_39	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	39	116878.70559999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_40	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	40	122087.95859999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_41	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	41	127405.74209999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_42	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	42	132832.0561	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_43	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	43	138366.9001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_44	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	44	144010.2741	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_45	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	45	149762.1776	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_46	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	46	155622.61109999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_47	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	47	161591.5751	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_48	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	48	167669.0691	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_49	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	49	173855.0931	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_50	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	50	180149.6476	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_51	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	51	186552.7326	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_52	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	52	193064.34759999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_53	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	53	199684.49259999997	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_54	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	54	206413.16709999996	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_55	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	55	213250.37159999995	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_56	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	56	220196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_57	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	57	227196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_58	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	58	234196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_59	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	59	241196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_60	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	60	248196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_61	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	61	255196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_62	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	62	262196.10659999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_101	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	101	0	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_102	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	102	2115.041	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_103	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	103	7221.6641500000005	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_104	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	104	13214.4211	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_105	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	105	19226.3636	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_106	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	106	25257.4918	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_107	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	107	31307.80555	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_108	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	108	37377.304950000005	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_109	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	109	43465.98990000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_110	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	110	49573.860400000005	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_111	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	111	55700.916450000004	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_112	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	112	61847.158050000005	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_113	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	113	68012.5853	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_114	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	114	74197.19810000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_115	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	115	80400.99665000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_116	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	116	86623.98075000002	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_117	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	117	92866.15035000001	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_118	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	118	99127.5056	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_119	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	119	105408.0464	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_120	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	120	111707.7727	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_121	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	121	118026.6847	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_122	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	122	124364.7825	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_123	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	123	130722.0658	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_124	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	124	137098.53465	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_125	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	125	143494.18899999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_126	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	126	149909.02899999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_127	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	127	156343.05449999997	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_128	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	128	162796.26549999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_129	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	129	169268.66249999998	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_130	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	130	175760.24499999997	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_131	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	131	182271.01249999995	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_132	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	132	188800.96599999996	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_133	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	133	195350.10549999995	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_134	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	134	201918.42999999993	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_135	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	135	208505.93949999995	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_136	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	136	215112.63499999995	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_137	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	137	221738.51649999994	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_138	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	138	228383.58299999993	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_139	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	139	235047.83499999993	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_140	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	140	241731.27299999993	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_141	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	141	248433.89649999992	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_142	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	142	255155.7049999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_143	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	143	261896.6994999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_144	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	144	268656.8799999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_145	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	145	275436.2454999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_146	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	146	282234.7964999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_147	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	147	289052.5334999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_148	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	148	295889.4559999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_149	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	149	302745.5634999999	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_150	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	150	309620.85699999984	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_151	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	151	316515.33649999986	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_152	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	152	323429.00149999984	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_153	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	153	330361.85199999984	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_154	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	154	337313.88799999986	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_155	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	155	344285.10949999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_156	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	156	351275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_157	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	157	358275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_158	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	158	365275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_159	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	159	372275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_160	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	160	379275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_161	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	161	386275.51649999985	calc_mag_energy_density	3
cp_mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density_162	mod_calc_mag_energy_density_Ni2MnGa_sample_1_mag_energy_density	mag_energy_density	Ni2MnGa_sample_1	162	393275.51649999985	calc_mag_energy_density	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_10_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_10_relative_permeability	relative_permeability	Ni2MnGa_sample_1	10	2.7751742798026418	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_102_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_102_relative_permeability	relative_permeability	Ni2MnGa_sample_1	102	34.66291580455197	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_103_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_103_relative_permeability	relative_permeability	Ni2MnGa_sample_1	103	24.806956469839257	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_104_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_104_relative_permeability	relative_permeability	Ni2MnGa_sample_1	104	16.9221964029922	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_105_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_105_relative_permeability	relative_permeability	Ni2MnGa_sample_1	105	12.97981736431641	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_106_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_106_relative_permeability	relative_permeability	Ni2MnGa_sample_1	106	10.614389622791661	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_107_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_107_relative_permeability	relative_permeability	Ni2MnGa_sample_1	107	9.037437662475462	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_108_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_108_relative_permeability	relative_permeability	Ni2MnGa_sample_1	108	7.911043746163116	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_109_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_109_relative_permeability	relative_permeability	Ni2MnGa_sample_1	109	7.066247811554989	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_11_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_11_relative_permeability	relative_permeability	Ni2MnGa_sample_1	11	2.6840247493235716	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_110_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_110_relative_permeability	relative_permeability	Ni2MnGa_sample_1	110	6.409184748969883	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_111_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_111_relative_permeability	relative_permeability	Ni2MnGa_sample_1	111	5.883533901002706	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_112_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_112_relative_permeability	relative_permeability	Ni2MnGa_sample_1	112	5.453456296210554	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_113_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_113_relative_permeability	relative_permeability	Ni2MnGa_sample_1	113	5.095058093267547	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_114_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_114_relative_permeability	relative_permeability	Ni2MnGa_sample_1	114	4.791798136607941	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_115_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_115_relative_permeability	relative_permeability	Ni2MnGa_sample_1	115	4.531861201427889	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_116_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_116_relative_permeability	relative_permeability	Ni2MnGa_sample_1	116	4.306582152899359	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_117_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_117_relative_permeability	relative_permeability	Ni2MnGa_sample_1	117	4.1094632838612135	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_118_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_118_relative_permeability	relative_permeability	Ni2MnGa_sample_1	118	3.9355347295690515	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_119_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_119_relative_permeability	relative_permeability	Ni2MnGa_sample_1	119	3.780931614409252	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_12_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_12_relative_permeability	relative_permeability	Ni2MnGa_sample_1	12	2.6094480777856557	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_120_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_120_relative_permeability	relative_permeability	Ni2MnGa_sample_1	120	3.6426024276033075	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_121_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_121_relative_permeability	relative_permeability	Ni2MnGa_sample_1	121	3.5181063982174123	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_122_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_122_relative_permeability	relative_permeability	Ni2MnGa_sample_1	122	3.4054669819543286	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_123_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_123_relative_permeability	relative_permeability	Ni2MnGa_sample_1	123	3.303067440278964	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_124_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_124_relative_permeability	relative_permeability	Ni2MnGa_sample_1	124	3.2095723103751324	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_125_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_125_relative_permeability	relative_permeability	Ni2MnGa_sample_1	125	3.1238683086635897	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_126_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_126_relative_permeability	relative_permeability	Ni2MnGa_sample_1	126	3.0450208499124622	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_127_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_127_relative_permeability	relative_permeability	Ni2MnGa_sample_1	127	2.9722382742198117	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_128_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_128_relative_permeability	relative_permeability	Ni2MnGa_sample_1	128	2.9048472951703888	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_129_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_129_relative_permeability	relative_permeability	Ni2MnGa_sample_1	129	2.8422699574816397	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_13_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_13_relative_permeability	relative_permeability	Ni2MnGa_sample_1	13	2.547300917820574	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_130_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_130_relative_permeability	relative_permeability	Ni2MnGa_sample_1	130	2.784008023841016	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_131_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_131_relative_permeability	relative_permeability	Ni2MnGa_sample_1	131	2.729630219109767	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_132_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_132_relative_permeability	relative_permeability	Ni2MnGa_sample_1	132	2.6787611732632346	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_133_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_133_relative_permeability	relative_permeability	Ni2MnGa_sample_1	133	2.6310709454082444	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_134_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_134_relative_permeability	relative_permeability	Ni2MnGa_sample_1	134	2.5862710343929507	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_135_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_135_relative_permeability	relative_permeability	Ni2MnGa_sample_1	135	2.54410641226091	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_136_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_136_relative_permeability	relative_permeability	Ni2MnGa_sample_1	136	2.504351651849663	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_137_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_137_relative_permeability	relative_permeability	Ni2MnGa_sample_1	137	2.466805047128937	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_138_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_138_relative_permeability	relative_permeability	Ni2MnGa_sample_1	138	2.4312879886093324	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_139_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_139_relative_permeability	relative_permeability	Ni2MnGa_sample_1	139	2.397640458379755	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_14_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_14_relative_permeability	relative_permeability	Ni2MnGa_sample_1	14	2.494714675742847	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_140_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_140_relative_permeability	relative_permeability	Ni2MnGa_sample_1	140	2.3657184425209254	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_141_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_141_relative_permeability	relative_permeability	Ni2MnGa_sample_1	141	2.335392328505491	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_142_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_142_relative_permeability	relative_permeability	Ni2MnGa_sample_1	142	2.306545537124956	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_143_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_143_relative_permeability	relative_permeability	Ni2MnGa_sample_1	143	2.2790727814283445	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_144_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_144_relative_permeability	relative_permeability	Ni2MnGa_sample_1	144	2.2528774581835815	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_145_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_145_relative_permeability	relative_permeability	Ni2MnGa_sample_1	145	2.227872831449944	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_146_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_146_relative_permeability	relative_permeability	Ni2MnGa_sample_1	146	2.2039796983040656	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_147_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_147_relative_permeability	relative_permeability	Ni2MnGa_sample_1	147	2.181125397034095	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_148_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_148_relative_permeability	relative_permeability	Ni2MnGa_sample_1	148	2.159243449903658	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_149_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_149_relative_permeability	relative_permeability	Ni2MnGa_sample_1	149	2.1382732505703217	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_15_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_15_relative_permeability	relative_permeability	Ni2MnGa_sample_1	15	2.449640753961938	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_150_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_150_relative_permeability	relative_permeability	Ni2MnGa_sample_1	150	2.1181593025559735	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_151_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_151_relative_permeability	relative_permeability	Ni2MnGa_sample_1	151	2.0988495941429255	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_152_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_152_relative_permeability	relative_permeability	Ni2MnGa_sample_1	152	2.080297285235915	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_153_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_153_relative_permeability	relative_permeability	Ni2MnGa_sample_1	153	2.0624583736333695	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_154_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_154_relative_permeability	relative_permeability	Ni2MnGa_sample_1	154	2.045292778656993	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_155_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_155_relative_permeability	relative_permeability	Ni2MnGa_sample_1	155	2.0287627990874846	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_156_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_156_relative_permeability	relative_permeability	Ni2MnGa_sample_1	156	2.0128340543747196	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_157_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_157_relative_permeability	relative_permeability	Ni2MnGa_sample_1	157	1.9947477319751712	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_158_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_158_relative_permeability	relative_permeability	Ni2MnGa_sample_1	158	1.9772960173791156	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_159_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_159_relative_permeability	relative_permeability	Ni2MnGa_sample_1	159	1.9604460860449928	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_16_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_16_relative_permeability	relative_permeability	Ni2MnGa_sample_1	16	2.410576688418484	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_160_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_160_relative_permeability	relative_permeability	Ni2MnGa_sample_1	160	1.9441673388238914	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_161_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_161_relative_permeability	relative_permeability	Ni2MnGa_sample_1	161	1.9284312165101598	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_162_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_162_relative_permeability	relative_permeability	Ni2MnGa_sample_1	162	1.9132110326329441	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_17_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_17_relative_permeability	relative_permeability	Ni2MnGa_sample_1	17	2.376395780280121	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_18_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_18_relative_permeability	relative_permeability	Ni2MnGa_sample_1	18	2.346236108640496	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_19_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_19_relative_permeability	relative_permeability	Ni2MnGa_sample_1	19	2.3194272021504236	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_2_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_2_relative_permeability	relative_permeability	Ni2MnGa_sample_1	2	8.979273436256566	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_20_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_20_relative_permeability	relative_permeability	Ni2MnGa_sample_1	20	2.295440662774236	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_21_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_21_relative_permeability	relative_permeability	Ni2MnGa_sample_1	21	2.2738526977558493	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_22_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_22_relative_permeability	relative_permeability	Ni2MnGa_sample_1	22	2.2543208051961834	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_23_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_23_relative_permeability	relative_permeability	Ni2MnGa_sample_1	23	2.2365643945422713	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_24_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_24_relative_permeability	relative_permeability	Ni2MnGa_sample_1	24	2.2203519157976324	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_25_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_25_relative_permeability	relative_permeability	Ni2MnGa_sample_1	25	2.205490675897926	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_26_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_26_relative_permeability	relative_permeability	Ni2MnGa_sample_1	26	2.1918182396944137	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_27_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_27_relative_permeability	relative_permeability	Ni2MnGa_sample_1	27	2.179197621175578	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_28_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_28_relative_permeability	relative_permeability	Ni2MnGa_sample_1	28	2.167511863287767	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_29_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_29_relative_permeability	relative_permeability	Ni2MnGa_sample_1	29	2.156660546599668	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_3_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_3_relative_permeability	relative_permeability	Ni2MnGa_sample_1	3	5.965402275982811	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_30_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_30_relative_permeability	relative_permeability	Ni2MnGa_sample_1	30	2.1465578709928814	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_31_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_31_relative_permeability	relative_permeability	Ni2MnGa_sample_1	31	2.1371286540400023	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_32_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_32_relative_permeability	relative_permeability	Ni2MnGa_sample_1	32	2.128307825006546	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_33_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_33_relative_permeability	relative_permeability	Ni2MnGa_sample_1	33	2.1200380739694413	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_34_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_34_relative_permeability	relative_permeability	Ni2MnGa_sample_1	34	2.1122697611158534	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_35_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_35_relative_permeability	relative_permeability	Ni2MnGa_sample_1	35	2.104958290812744	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_36_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_36_relative_permeability	relative_permeability	Ni2MnGa_sample_1	36	2.098064618812669	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_37_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_37_relative_permeability	relative_permeability	Ni2MnGa_sample_1	37	2.0915540170123967	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_38_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_38_relative_permeability	relative_permeability	Ni2MnGa_sample_1	38	2.0853951675692883	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_39_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_39_relative_permeability	relative_permeability	Ni2MnGa_sample_1	39	2.0795605518650992	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_4_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_4_relative_permeability	relative_permeability	Ni2MnGa_sample_1	4	4.598161706191309	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_40_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_40_relative_permeability	relative_permeability	Ni2MnGa_sample_1	40	2.074025147222664	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_41_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_41_relative_permeability	relative_permeability	Ni2MnGa_sample_1	41	2.0687667117618975	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_42_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_42_relative_permeability	relative_permeability	Ni2MnGa_sample_1	42	2.0637645912508784	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_43_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_43_relative_permeability	relative_permeability	Ni2MnGa_sample_1	43	2.0590006669546694	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_44_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_44_relative_permeability	relative_permeability	Ni2MnGa_sample_1	44	2.054458320532704	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_45_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_45_relative_permeability	relative_permeability	Ni2MnGa_sample_1	45	2.0501222635394214	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_46_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_46_relative_permeability	relative_permeability	Ni2MnGa_sample_1	46	2.0459792738783666	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_47_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_47_relative_permeability	relative_permeability	Ni2MnGa_sample_1	47	2.0420162412029708	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_48_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_48_relative_permeability	relative_permeability	Ni2MnGa_sample_1	48	2.0382218482158883	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_49_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_49_relative_permeability	relative_permeability	Ni2MnGa_sample_1	49	2.0345855549366014	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_5_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_5_relative_permeability	relative_permeability	Ni2MnGa_sample_1	5	3.9145412223460134	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_50_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_50_relative_permeability	relative_permeability	Ni2MnGa_sample_1	50	2.031097844198956	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_51_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_51_relative_permeability	relative_permeability	Ni2MnGa_sample_1	51	2.0277494827311795	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_52_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_52_relative_permeability	relative_permeability	Ni2MnGa_sample_1	52	2.024532429556257	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_53_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_53_relative_permeability	relative_permeability	Ni2MnGa_sample_1	53	2.0214391091957538	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_54_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_54_relative_permeability	relative_permeability	Ni2MnGa_sample_1	54	2.018462367754858	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_55_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_55_relative_permeability	relative_permeability	Ni2MnGa_sample_1	55	2.015596170737027	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_56_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_56_relative_permeability	relative_permeability	Ni2MnGa_sample_1	56	2.0128340543747196	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_57_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_57_relative_permeability	relative_permeability	Ni2MnGa_sample_1	57	1.9947477319751712	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_58_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_58_relative_permeability	relative_permeability	Ni2MnGa_sample_1	58	1.9772960173791156	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_59_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_59_relative_permeability	relative_permeability	Ni2MnGa_sample_1	59	1.9604460860449928	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_6_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_6_relative_permeability	relative_permeability	Ni2MnGa_sample_1	6	3.504369091198472	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_60_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_60_relative_permeability	relative_permeability	Ni2MnGa_sample_1	60	1.9441673388238914	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_61_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_61_relative_permeability	relative_permeability	Ni2MnGa_sample_1	61	1.9284312165101598	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_62_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_62_relative_permeability	relative_permeability	Ni2MnGa_sample_1	62	1.9132110326329441	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_7_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_7_relative_permeability	relative_permeability	Ni2MnGa_sample_1	7	3.2309210037667784	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_8_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_8_relative_permeability	relative_permeability	Ni2MnGa_sample_1	8	3.0356010550010235	calc_rel_perm	3
cp_mod_relative_permeability_calc_rel_perm_Ni2MnGa_sample_1_9_relative_permeability	mod_calc_rel_perm_relative_permeability_Ni2MnGa_sample_1_9_relative_permeability	relative_permeability	Ni2MnGa_sample_1	9	2.889110894477161	calc_rel_perm	3
cp_mod_A01_matrix_model_1_M420_0_1	mod_matrix_model_1_A01_M420_0	A01	M420	0	1.35929071606055e-11	matrix_model_1	3
cp_mod_A01_matrix_model_1_M420_0_2	mod_matrix_model_1_A01_M420_0	A01	M420	0	-7.50709283939466e-12	matrix_model_1	3
cp_mod_A01_matrix_model_1_M420_0_3	mod_matrix_model_1_A01_M420_0	A01	M420	0	-2.49051276259319e-12	matrix_model_1	3
cp_mod_A01_matrix_model_1_M420_0_4	mod_matrix_model_1_A01_M420_0	A01	M420	0	9.80395019200366e-12	matrix_model_1	3
cp_mod_A01_matrix_model_1_M420_0_5	mod_matrix_model_1_A01_M420_0	A01	M420	0	2.55437514117913e-11	matrix_model_1	3
cp_mod_B01_matrix_model_2_M420_0_1	mod_matrix_model_2_B01_M420_0	B01	M420	0	-0.011	matrix_model_2	3
cp_mod_B01_matrix_model_2_M420_0_2	mod_matrix_model_2_B01_M420_0	B01	M420	0	0.0250000000000002	matrix_model_2	3
cp_mod_B01_matrix_model_2_M420_0_3	mod_matrix_model_2_B01_M420_0	B01	M420	0	0.0370595211203975	matrix_model_2	3
cp_mod_C01_matrix_model_3_M420_0_1	mod_matrix_model_3_C01_M420_0	C01	M420	0	70589564.0388528	matrix_model_3	3
cp_mod_C01_matrix_model_3_M420_0_2	mod_matrix_model_3_C01_M420_0	C01	M420	0	72002066.8080778	matrix_model_3	3
cp_mod_A10_matrix_model_4_M420_0_1	mod_matrix_model_4_A10_M420_0	A10	M420	0	120993245639.767	matrix_model_4	3
cp_mod_A10_matrix_model_4_M420_0_2	mod_matrix_model_4_A10_M420_0	A10	M420	0	73599880710.8577	matrix_model_4	3
cp_mod_A10_matrix_model_4_M420_0_3	mod_matrix_model_4_A10_M420_0	A10	M420	0	66426959648.0118	matrix_model_4	3
cp_mod_A10_matrix_model_4_M420_0_4	mod_matrix_model_4_A10_M420_0	A10	M420	0	99451550706.5426	matrix_model_4	3
cp_mod_A10_matrix_model_4_M420_0_5	mod_matrix_model_4_A10_M420_0	A10	M420	0	22222222222.2223	matrix_model_4	3
cp_mod_B10_matrix_model_5_M420_0_1	mod_matrix_model_5_B10_M420_0	B10	M420	0	-6.66439756424805	matrix_model_5	3
cp_mod_B10_matrix_model_5_M420_0_2	mod_matrix_model_5_B10_M420_0	B10	M420	0	14.234253276912	matrix_model_5	3
cp_mod_B10_matrix_model_5_M420_0_3	mod_matrix_model_5_B10_M420_0	B10	M420	0	11.6666666666666	matrix_model_5	3
cp_mod_C10_matrix_model_6_M420_0_1	mod_matrix_model_6_C10_M420_0	C10	M420	0	8.04139999999999e-09	matrix_model_6	3
cp_mod_C10_matrix_model_6_M420_0_2	mod_matrix_model_6_C10_M420_0	C10	M420	0	6.90989778098874e-09	matrix_model_6	3
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_610000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1224	610000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_610000	\N	magnetic_stress	Ni2MnGa_sample_1	1224	1933421.2960250007	calc_magnetic_stress	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_10_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_10_flux_density	flux_density	Ni2MnGa_sample_1	10	0.31385556	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_102_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_102_flux_density	flux_density	Ni2MnGa_sample_1	102	0.4355742000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_103_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_103_flux_density	flux_density	Ni2MnGa_sample_1	103	0.6234484300000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_104_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_104_flux_density	flux_density	Ni2MnGa_sample_1	104	0.6379329599999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_105_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_105_flux_density	flux_density	Ni2MnGa_sample_1	105	0.6524175400000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_106_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_106_flux_density	flux_density	Ni2MnGa_sample_1	106	0.6669021	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_107_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_107_flux_density	flux_density	Ni2MnGa_sample_1	107	0.68138665	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_108_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_108_flux_density	flux_density	Ni2MnGa_sample_1	108	0.69587123	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_109_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_109_flux_density	flux_density	Ni2MnGa_sample_1	109	0.7103557599999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_11_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_11_flux_density	flux_density	Ni2MnGa_sample_1	11	0.33727455	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_110_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_110_flux_density	flux_density	Ni2MnGa_sample_1	110	0.7248403399999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_111_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_111_flux_density	flux_density	Ni2MnGa_sample_1	111	0.73932487	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_112_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_112_flux_density	flux_density	Ni2MnGa_sample_1	112	0.7538094500000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_113_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_113_flux_density	flux_density	Ni2MnGa_sample_1	113	0.768294	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_114_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_114_flux_density	flux_density	Ni2MnGa_sample_1	114	0.7827785599999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_115_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_115_flux_density	flux_density	Ni2MnGa_sample_1	115	0.7972631499999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_116_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_116_flux_density	flux_density	Ni2MnGa_sample_1	116	0.81174767	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_117_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_117_flux_density	flux_density	Ni2MnGa_sample_1	117	0.8262322500000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_118_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_118_flux_density	flux_density	Ni2MnGa_sample_1	118	0.8407167999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_119_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_119_flux_density	flux_density	Ni2MnGa_sample_1	119	0.8552013599999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_12_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_12_flux_density	flux_density	Ni2MnGa_sample_1	12	0.3606935700000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_120_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_120_flux_density	flux_density	Ni2MnGa_sample_1	120	0.8696859	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_121_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_121_flux_density	flux_density	Ni2MnGa_sample_1	121	0.8841705000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_122_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_122_flux_density	flux_density	Ni2MnGa_sample_1	122	0.8986550600000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_123_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_123_flux_density	flux_density	Ni2MnGa_sample_1	123	0.9131396	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_124_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_124_flux_density	flux_density	Ni2MnGa_sample_1	124	0.92762417	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_125_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_125_flux_density	flux_density	Ni2MnGa_sample_1	125	0.9421086999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_126_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_126_flux_density	flux_density	Ni2MnGa_sample_1	126	0.9565932999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_127_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_127_flux_density	flux_density	Ni2MnGa_sample_1	127	0.9710778	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_128_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_128_flux_density	flux_density	Ni2MnGa_sample_1	128	0.9855623999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_129_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_129_flux_density	flux_density	Ni2MnGa_sample_1	129	1.000047	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_13_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_13_flux_density	flux_density	Ni2MnGa_sample_1	13	0.3841126	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_130_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_130_flux_density	flux_density	Ni2MnGa_sample_1	130	1.0145315	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_131_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_131_flux_density	flux_density	Ni2MnGa_sample_1	131	1.029016	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_132_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_132_flux_density	flux_density	Ni2MnGa_sample_1	132	1.0435007	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_133_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_133_flux_density	flux_density	Ni2MnGa_sample_1	133	1.0579851999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_134_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_134_flux_density	flux_density	Ni2MnGa_sample_1	134	1.0724696999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_135_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_135_flux_density	flux_density	Ni2MnGa_sample_1	135	1.0869542000000003	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_136_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_136_flux_density	flux_density	Ni2MnGa_sample_1	136	1.1014389000000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_137_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_137_flux_density	flux_density	Ni2MnGa_sample_1	137	1.1159234	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_138_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_138_flux_density	flux_density	Ni2MnGa_sample_1	138	1.1304079000000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_139_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_139_flux_density	flux_density	Ni2MnGa_sample_1	139	1.1448924999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_14_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_14_flux_density	flux_density	Ni2MnGa_sample_1	14	0.4075316	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_140_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_140_flux_density	flux_density	Ni2MnGa_sample_1	140	1.1593771	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_141_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_141_flux_density	flux_density	Ni2MnGa_sample_1	141	1.1738616	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_142_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_142_flux_density	flux_density	Ni2MnGa_sample_1	142	1.1883461	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_143_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_143_flux_density	flux_density	Ni2MnGa_sample_1	143	1.2028308	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_144_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_144_flux_density	flux_density	Ni2MnGa_sample_1	144	1.2173153	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_145_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_145_flux_density	flux_density	Ni2MnGa_sample_1	145	1.2317998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_146_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_146_flux_density	flux_density	Ni2MnGa_sample_1	146	1.2462844	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_147_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_147_flux_density	flux_density	Ni2MnGa_sample_1	147	1.260769	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_148_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_148_flux_density	flux_density	Ni2MnGa_sample_1	148	1.2752535000000003	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_149_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_149_flux_density	flux_density	Ni2MnGa_sample_1	149	1.2897379999999996	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_15_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_15_flux_density	flux_density	Ni2MnGa_sample_1	15	0.4309505999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_150_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_150_flux_density	flux_density	Ni2MnGa_sample_1	150	1.3042226999999997	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_151_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_151_flux_density	flux_density	Ni2MnGa_sample_1	151	1.3187072	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_152_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_152_flux_density	flux_density	Ni2MnGa_sample_1	152	1.3331918	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_153_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_153_flux_density	flux_density	Ni2MnGa_sample_1	153	1.3476762999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_154_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_154_flux_density	flux_density	Ni2MnGa_sample_1	154	1.3621609	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_155_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_155_flux_density	flux_density	Ni2MnGa_sample_1	155	1.3766454	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_156_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_156_flux_density	flux_density	Ni2MnGa_sample_1	156	1.39113	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_157_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_157_flux_density	flux_density	Ni2MnGa_sample_1	157	1.4036959999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_158_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_158_flux_density	flux_density	Ni2MnGa_sample_1	158	1.4162620000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_159_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_159_flux_density	flux_density	Ni2MnGa_sample_1	159	1.428828	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_16_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_16_flux_density	flux_density	Ni2MnGa_sample_1	16	0.4543696	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_160_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_160_flux_density	flux_density	Ni2MnGa_sample_1	160	1.441394	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_161_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_161_flux_density	flux_density	Ni2MnGa_sample_1	161	1.4539600000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_162_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_162_flux_density	flux_density	Ni2MnGa_sample_1	162	1.466526	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_17_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_17_flux_density	flux_density	Ni2MnGa_sample_1	17	0.47778862999999994	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_18_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_18_flux_density	flux_density	Ni2MnGa_sample_1	18	0.50120765	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_19_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_19_flux_density	flux_density	Ni2MnGa_sample_1	19	0.5246266	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_2_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_2_flux_density	flux_density	Ni2MnGa_sample_1	2	0.11283355000000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_20_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_20_flux_density	flux_density	Ni2MnGa_sample_1	20	0.54804564	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_21_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_21_flux_density	flux_density	Ni2MnGa_sample_1	21	0.57146466	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_22_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_22_flux_density	flux_density	Ni2MnGa_sample_1	22	0.5948837	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_23_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_23_flux_density	flux_density	Ni2MnGa_sample_1	23	0.6183027	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_24_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_24_flux_density	flux_density	Ni2MnGa_sample_1	24	0.6417216700000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_25_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_25_flux_density	flux_density	Ni2MnGa_sample_1	25	0.6651407	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_26_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_26_flux_density	flux_density	Ni2MnGa_sample_1	26	0.6885597	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_27_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_27_flux_density	flux_density	Ni2MnGa_sample_1	27	0.7119787300000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_28_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_28_flux_density	flux_density	Ni2MnGa_sample_1	28	0.7353977600000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_29_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_29_flux_density	flux_density	Ni2MnGa_sample_1	29	0.7588167000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_3_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_3_flux_density	flux_density	Ni2MnGa_sample_1	3	0.14992249	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_30_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_30_flux_density	flux_density	Ni2MnGa_sample_1	30	0.78223574	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_31_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_31_flux_density	flux_density	Ni2MnGa_sample_1	31	0.80565476	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_32_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_32_flux_density	flux_density	Ni2MnGa_sample_1	32	0.8290738	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_33_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_33_flux_density	flux_density	Ni2MnGa_sample_1	33	0.8524927499999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_34_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_34_flux_density	flux_density	Ni2MnGa_sample_1	34	0.8759117999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_35_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_35_flux_density	flux_density	Ni2MnGa_sample_1	35	0.8993308	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_36_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_36_flux_density	flux_density	Ni2MnGa_sample_1	36	0.9227498	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_37_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_37_flux_density	flux_density	Ni2MnGa_sample_1	37	0.94616884	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_38_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_38_flux_density	flux_density	Ni2MnGa_sample_1	38	0.9695878000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_39_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_39_flux_density	flux_density	Ni2MnGa_sample_1	39	0.9930067999999997	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_4_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_4_flux_density	flux_density	Ni2MnGa_sample_1	4	0.17334149999999995	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_40_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_40_flux_density	flux_density	Ni2MnGa_sample_1	40	1.0164257999999997	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_41_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_41_flux_density	flux_density	Ni2MnGa_sample_1	41	1.0398449	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_42_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_42_flux_density	flux_density	Ni2MnGa_sample_1	42	1.0632639000000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_43_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_43_flux_density	flux_density	Ni2MnGa_sample_1	43	1.0866828999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_44_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_44_flux_density	flux_density	Ni2MnGa_sample_1	44	1.1101019	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_45_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_45_flux_density	flux_density	Ni2MnGa_sample_1	45	1.1335208	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_46_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_46_flux_density	flux_density	Ni2MnGa_sample_1	46	1.1569399	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_47_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_47_flux_density	flux_density	Ni2MnGa_sample_1	47	1.1803589000000003	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_48_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_48_flux_density	flux_density	Ni2MnGa_sample_1	48	1.2037779000000002	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_49_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_49_flux_density	flux_density	Ni2MnGa_sample_1	49	1.2271969	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_5_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_5_flux_density	flux_density	Ni2MnGa_sample_1	5	0.19676050000000003	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_50_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_50_flux_density	flux_density	Ni2MnGa_sample_1	50	1.250616	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_51_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_51_flux_density	flux_density	Ni2MnGa_sample_1	51	1.274035	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_52_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_52_flux_density	flux_density	Ni2MnGa_sample_1	52	1.297454	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_53_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_53_flux_density	flux_density	Ni2MnGa_sample_1	53	1.3208729999999997	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_54_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_54_flux_density	flux_density	Ni2MnGa_sample_1	54	1.3442919	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_55_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_55_flux_density	flux_density	Ni2MnGa_sample_1	55	1.3677110000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_56_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_56_flux_density	flux_density	Ni2MnGa_sample_1	56	1.39113	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_57_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_57_flux_density	flux_density	Ni2MnGa_sample_1	57	1.4036959999999998	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_58_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_58_flux_density	flux_density	Ni2MnGa_sample_1	58	1.4162620000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_59_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_59_flux_density	flux_density	Ni2MnGa_sample_1	59	1.428828	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_6_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_6_flux_density	flux_density	Ni2MnGa_sample_1	6	0.22017951	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_60_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_60_flux_density	flux_density	Ni2MnGa_sample_1	60	1.441394	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_61_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_61_flux_density	flux_density	Ni2MnGa_sample_1	61	1.4539600000000001	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_62_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_62_flux_density	flux_density	Ni2MnGa_sample_1	62	1.466526	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_7_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_7_flux_density	flux_density	Ni2MnGa_sample_1	7	0.24359851999999999	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_8_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_8_flux_density	flux_density	Ni2MnGa_sample_1	8	0.26701754	calc_flux_density	4
cp_mod_flux_density_calc_flux_density_Ni2MnGa_sample_1_9_flux_density	mod_calc_flux_density_flux_density_Ni2MnGa_sample_1_9_flux_density	flux_density	Ni2MnGa_sample_1	9	0.29043654	calc_flux_density	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_160000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1134	160000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_160000	\N	magnetic_stress	Ni2MnGa_sample_1	1134	931210.3472750011	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_130000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1128	130000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_130000	\N	magnetic_stress	Ni2MnGa_sample_1	1128	771093.8005500009	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_290000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1160	290000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_290000	\N	magnetic_stress	Ni2MnGa_sample_1	1160	1487994.242600001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_10000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1104	10000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_10000	\N	magnetic_stress	Ni2MnGa_sample_1	1104	23802.122937500026	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_200000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1142	200000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_200000	\N	magnetic_stress	Ni2MnGa_sample_1	1142	1126249.4350250012	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_0magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1102	0	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_0	\N	magnetic_stress	Ni2MnGa_sample_1	1102	0	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_80000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1118	80000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_80000	\N	magnetic_stress	Ni2MnGa_sample_1	1118	477876.2687250006	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_590000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1220	590000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_590000	\N	magnetic_stress	Ni2MnGa_sample_1	1220	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_500000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1202	500000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_500000	\N	magnetic_stress	Ni2MnGa_sample_1	1202	1916948.4075250002	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_360000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1174	360000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_360000	\N	magnetic_stress	Ni2MnGa_sample_1	1174	1695552.724725001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_150000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1132	150000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_150000	\N	magnetic_stress	Ni2MnGa_sample_1	1132	879155.9976375011	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_420000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1186	420000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_420000	\N	magnetic_stress	Ni2MnGa_sample_1	1186	1822064.5411500004	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_220000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1146	220000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_220000	\N	magnetic_stress	Ni2MnGa_sample_1	1146	1215861.9877000013	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_580000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1218	580000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_580000	\N	magnetic_stress	Ni2MnGa_sample_1	1218	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_270000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1156	270000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_270000	\N	magnetic_stress	Ni2MnGa_sample_1	1156	1416831.318850001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_110000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1124	110000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_110000	\N	magnetic_stress	Ni2MnGa_sample_1	1124	657760.2846000007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_70000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1116	70000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_70000	\N	magnetic_stress	Ni2MnGa_sample_1	1116	415279.2666125005	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_530000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1208	530000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_530000	\N	magnetic_stress	Ni2MnGa_sample_1	1208	1930785.6332750004	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_370000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1176	370000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_370000	\N	magnetic_stress	Ni2MnGa_sample_1	1176	1719932.6059000008	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_460000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1194	460000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_460000	\N	magnetic_stress	Ni2MnGa_sample_1	1194	1880049.1364000007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_190000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1140	190000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_190000	\N	magnetic_stress	Ni2MnGa_sample_1	1140	1079466.410150001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_20000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1106	20000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_20000	\N	magnetic_stress	Ni2MnGa_sample_1	1106	82526.78395000008	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_250000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1152	250000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_250000	\N	magnetic_stress	Ni2MnGa_sample_1	1152	1340397.085825001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_100000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1122	100000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_100000	\N	magnetic_stress	Ni2MnGa_sample_1	1122	599116.7773500007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_560000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1214	560000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_560000	\N	magnetic_stress	Ni2MnGa_sample_1	1214	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_520000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1206	520000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_520000	\N	magnetic_stress	Ni2MnGa_sample_1	1206	1927491.05115	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_60000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1114	60000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_60000	\N	magnetic_stress	Ni2MnGa_sample_1	1114	351364.4323875004	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_350000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1172	350000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_350000	\N	magnetic_stress	Ni2MnGa_sample_1	1172	1669855.004800001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_570000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1216	570000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_570000	\N	magnetic_stress	Ni2MnGa_sample_1	1216	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_280000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1158	280000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_280000	\N	magnetic_stress	Ni2MnGa_sample_1	1158	1453071.695675001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_380000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1178	380000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_380000	\N	magnetic_stress	Ni2MnGa_sample_1	1178	1742994.658650001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_240000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1150	240000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_240000	\N	magnetic_stress	Ni2MnGa_sample_1	1150	1300203.217825001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_450000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1192	450000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_450000	\N	magnetic_stress	Ni2MnGa_sample_1	1192	1867529.7346500005	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_540000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1210	540000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_540000	\N	magnetic_stress	Ni2MnGa_sample_1	1210	1932762.3840250003	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_180000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1138	180000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_180000	\N	magnetic_stress	Ni2MnGa_sample_1	1138	1031365.552425001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_300000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1162	300000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_300000	\N	magnetic_stress	Ni2MnGa_sample_1	1162	1521598.9463500008	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_50000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1112	50000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_50000	\N	magnetic_stress	Ni2MnGa_sample_1	1112	286131.7667875003	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_340000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1170	340000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_340000	\N	magnetic_stress	Ni2MnGa_sample_1	1170	1642839.457925001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_170000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1136	170000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_170000	\N	magnetic_stress	Ni2MnGa_sample_1	1136	981946.864062501	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_40000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1110	40000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_40000	\N	magnetic_stress	Ni2MnGa_sample_1	1110	219581.26907500022	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_430000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1188	430000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_430000	\N	magnetic_stress	Ni2MnGa_sample_1	1188	1838537.437025	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_210000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1144	210000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_210000	\N	magnetic_stress	Ni2MnGa_sample_1	1144	1171714.6285250012	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_470000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1196	470000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_470000	\N	magnetic_stress	Ni2MnGa_sample_1	1196	1891250.7067750003	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_310000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1164	310000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_310000	\N	magnetic_stress	Ni2MnGa_sample_1	1164	1553885.8216750007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_440000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1190	440000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_440000	\N	magnetic_stress	Ni2MnGa_sample_1	1190	1853692.5015250004	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_390000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1180	390000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_390000	\N	magnetic_stress	Ni2MnGa_sample_1	1180	1764738.887400001	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_230000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1148	230000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_230000	\N	magnetic_stress	Ni2MnGa_sample_1	1148	1258691.5184500013	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_600000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1222	600000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_600000	\N	magnetic_stress	Ni2MnGa_sample_1	1222	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_120000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1126	120000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_120000	\N	magnetic_stress	Ni2MnGa_sample_1	1126	715085.9590000008	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_30000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1108	30000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_30000	\N	magnetic_stress	Ni2MnGa_sample_1	1108	151712.94220000017	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_490000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1200	490000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_490000	\N	magnetic_stress	Ni2MnGa_sample_1	1200	1909700.3386499998	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_410000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1184	410000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_410000	\N	magnetic_stress	Ni2MnGa_sample_1	1184	1804273.8212750005	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_320000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1166	320000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_320000	\N	magnetic_stress	Ni2MnGa_sample_1	1166	1584854.8737375007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_90000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1120	90000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_90000	\N	magnetic_stress	Ni2MnGa_sample_1	1120	539155.4387250006	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_480000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1198	480000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_480000	\N	magnetic_stress	Ni2MnGa_sample_1	1198	1901134.4384	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_550000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1212	550000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_550000	\N	magnetic_stress	Ni2MnGa_sample_1	1212	1933421.2960250007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_510000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1204	510000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_510000	\N	magnetic_stress	Ni2MnGa_sample_1	1204	1922878.645025	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_140000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1130	140000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_140000	\N	magnetic_stress	Ni2MnGa_sample_1	1130	825783.814412501	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_330000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1168	330000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_330000	\N	magnetic_stress	Ni2MnGa_sample_1	1168	1614506.0870500007	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_260000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1154	260000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_260000	\N	magnetic_stress	Ni2MnGa_sample_1	1154	1379273.1202375009	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_400000magnetic_field_strength	\N	magnetic_field_strength	Ni2MnGa_sample_1	1182	400000	calc_magnetic_stress	4
cp_mod__magnetic_stress_calc_magnetic_stress_Ni2MnGa_sample_1_400000	\N	magnetic_stress	Ni2MnGa_sample_1	1182	1785165.2774000007	calc_magnetic_stress	4
cp_mod_calc_max_magn_stress_max_magnetic_stress_Ni2MnGa_sample_1_0	mod_calc_max_magn_stress_max_magnetic_stress_Ni2MnGa_sample_1	max_magnetic_stress	Ni2MnGa_sample_1	0	1933421.2960250007	calc_max_magn_stress	5
cp_mod_maximal_blocking_stress_load_calc_max_block_load_hold_stick_Ni2MnGa_sample_1_0_maximal_blocking_stress_load	mod_calc_max_block_load_hold_maximal_blocking_stress_load_stick_Ni2MnGa_sample_1_0_maximal_blocking_stress_load	maximal_blocking_stress_load	stick_Ni2MnGa_sample_1	0	1633421.2960250007	calc_max_block_load_hold	6
cp_mod_maximal_blocking_stress_load_calc_max_block_load_hold_stick_Ni2MnGa_sample_1_0_maximal_blocking_stress_hold	mod_calc_max_block_load_hold_maximal_blocking_stress_load_stick_Ni2MnGa_sample_1_0_maximal_blocking_stress_hold	maximal_blocking_stress_hold	stick_Ni2MnGa_sample_1	0	2233421.2960250005	calc_max_block_load_hold	6
\.


--
-- Name: characteristic_curve characteristic_curve_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.characteristic_curve
    ADD CONSTRAINT characteristic_curve_pkey PRIMARY KEY (curve_id);


--
-- Name: component_parts component_parts_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.component_parts
    ADD CONSTRAINT component_parts_pkey PRIMARY KEY (concr_component_id, has_part);


--
-- Name: component component_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.component
    ADD CONSTRAINT component_pkey PRIMARY KEY (component_id);


--
-- Name: concr_meas concr_meas_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_meas
    ADD CONSTRAINT concr_meas_pkey PRIMARY KEY (concr_meas_id);


--
-- Name: concr_model concr_model_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_model
    ADD CONSTRAINT concr_model_pkey PRIMARY KEY (concr_model_id);


--
-- Name: concr_param_matrix_additional concr_param_matrix_additional_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_param_matrix_additional
    ADD CONSTRAINT concr_param_matrix_additional_pkey PRIMARY KEY (concr_param_id, index_id);


--
-- Name: concr_param_matrix concr_param_matrix_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_param_matrix
    ADD CONSTRAINT concr_param_matrix_pkey PRIMARY KEY (concr_param_id, index_id);


--
-- Name: concr_param concr_param_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_param
    ADD CONSTRAINT concr_param_pkey PRIMARY KEY (concr_param_id);


--
-- Name: concr_component conr_component_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_component
    ADD CONSTRAINT conr_component_pkey PRIMARY KEY (concr_component_id);


--
-- Name: constants constants_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.constants
    ADD CONSTRAINT constants_pkey PRIMARY KEY (param_id, value);


--
-- Name: curve_condition curve_condition_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.curve_condition
    ADD CONSTRAINT curve_condition_pkey PRIMARY KEY (curve_id, param_id, in_number);


--
-- Name: curve_params curve_params_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.curve_params
    ADD CONSTRAINT curve_params_pkey PRIMARY KEY (curve_id, param_id);


--
-- Name: material_condition material_condition_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.material_condition
    ADD CONSTRAINT material_condition_pkey PRIMARY KEY (mod_id, mat_type);


--
-- Name: material material_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.material
    ADD CONSTRAINT material_pkey PRIMARY KEY (mat_id);


--
-- Name: material_types material_types_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.material_types
    ADD CONSTRAINT material_types_pkey PRIMARY KEY (mat_type);


--
-- Name: measinput measinput_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.measinput
    ADD CONSTRAINT measinput_pkey PRIMARY KEY (meas_id, param_id);


--
-- Name: measoutput measoutput_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.measoutput
    ADD CONSTRAINT measoutput_pkey PRIMARY KEY (meas_id, param_id);


--
-- Name: measurement measurement_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.measurement
    ADD CONSTRAINT measurement_pkey PRIMARY KEY (meas_id);


--
-- Name: model_condition_compare mod_con_compare_pk; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.model_condition_compare
    ADD CONSTRAINT mod_con_compare_pk PRIMARY KEY (mod_id, in_number, compare_position);


--
-- Name: model_condition model_condition_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.model_condition
    ADD CONSTRAINT model_condition_pkey PRIMARY KEY (mod_id, param_id);


--
-- Name: model model_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.model
    ADD CONSTRAINT model_pkey PRIMARY KEY (mod_id);


--
-- Name: modelinput_matrix modelinput_matrix_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modelinput_matrix
    ADD CONSTRAINT modelinput_matrix_pkey PRIMARY KEY (mod_id, param_id, index_id);


--
-- Name: modelinput modelinput_mod_id_in_number_key; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modelinput
    ADD CONSTRAINT modelinput_mod_id_in_number_key UNIQUE (mod_id, in_number);


--
-- Name: modelinput modelinput_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modelinput
    ADD CONSTRAINT modelinput_pkey PRIMARY KEY (mod_id, param_id);


--
-- Name: modeloutput_pek_bits modeloutput_pek_bits_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modeloutput_pek_bits
    ADD CONSTRAINT modeloutput_pek_bits_pkey PRIMARY KEY (mod_id);


--
-- Name: param_piezo param_piezo_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.param_piezo
    ADD CONSTRAINT param_piezo_pkey PRIMARY KEY (param_id, index_id);


--
-- Name: param param_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.param
    ADD CONSTRAINT param_pkey PRIMARY KEY (param_id);


--
-- Name: param_semantic param_semantic_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.param_semantic
    ADD CONSTRAINT param_semantic_pkey PRIMARY KEY (param_id, value);


--
-- Name: modeltype pk; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modeltype
    ADD CONSTRAINT pk PRIMARY KEY (modeltype);


--
-- Name: modeloutput pkk; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modeloutput
    ADD CONSTRAINT pkk PRIMARY KEY (mod_id, out_number);


--
-- Name: sample sample_pkey; Type: CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.sample
    ADD CONSTRAINT sample_pkey PRIMARY KEY (sample_id);


--
-- Name: concr_param_concr_param_id_idx; Type: INDEX; Schema: public; Owner: mena
--

CREATE INDEX concr_param_concr_param_id_idx ON public.concr_param USING btree (concr_param_id);


--
-- Name: fki_fks component parts; Type: INDEX; Schema: public; Owner: mena
--

CREATE INDEX "fki_fks component parts" ON public.component_parts USING btree (concr_component_id);


--
-- Name: ftk; Type: INDEX; Schema: public; Owner: mena
--

CREATE UNIQUE INDEX ftk ON public.concr_parameter USING btree (concr_param_id);


--
-- Name: param_semantic_unique; Type: INDEX; Schema: public; Owner: mena
--

CREATE UNIQUE INDEX param_semantic_unique ON public.param_semantic USING btree (param_id, value);


--
-- Name: concr_param_matrix concr_param_id; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_param_matrix
    ADD CONSTRAINT concr_param_id FOREIGN KEY (concr_param_id) REFERENCES public.concr_param(concr_param_id);


--
-- Name: concr_param concr_param_param_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_param
    ADD CONSTRAINT concr_param_param_id_fkey FOREIGN KEY (param_id) REFERENCES public.param(param_id) ON UPDATE CASCADE;


--
-- Name: constants constants_param_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.constants
    ADD CONSTRAINT constants_param_id_fkey FOREIGN KEY (param_id) REFERENCES public.param(param_id) ON UPDATE CASCADE;


--
-- Name: curve_params curve_id; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.curve_params
    ADD CONSTRAINT curve_id FOREIGN KEY (curve_id) REFERENCES public.characteristic_curve(curve_id);


--
-- Name: material fk_mat_type; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.material
    ADD CONSTRAINT fk_mat_type FOREIGN KEY (mat_type) REFERENCES public.material_types(mat_type) ON UPDATE CASCADE;


--
-- Name: component_parts fks component parts; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.component_parts
    ADD CONSTRAINT "fks component parts" FOREIGN KEY (concr_component_id) REFERENCES public.concr_component(concr_component_id) NOT VALID;


--
-- Name: concr_component ftk; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_component
    ADD CONSTRAINT ftk FOREIGN KEY (component_id) REFERENCES public.component(component_id);


--
-- Name: model_condition ftk; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.model_condition
    ADD CONSTRAINT ftk FOREIGN KEY (mod_id) REFERENCES public.model(mod_id);


--
-- Name: curve_condition ftk; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.curve_condition
    ADD CONSTRAINT ftk FOREIGN KEY (curve_id) REFERENCES public.characteristic_curve(curve_id);


--
-- Name: material_condition ftk; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.material_condition
    ADD CONSTRAINT ftk FOREIGN KEY (mod_id) REFERENCES public.model(mod_id);


--
-- Name: modelinput_matrix ftk; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modelinput_matrix
    ADD CONSTRAINT ftk FOREIGN KEY (mod_id, param_id) REFERENCES public.modelinput(mod_id, param_id);


--
-- Name: model_condition ftk1; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.model_condition
    ADD CONSTRAINT ftk1 FOREIGN KEY (param_id) REFERENCES public.param(param_id);


--
-- Name: curve_condition ftk1; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.curve_condition
    ADD CONSTRAINT ftk1 FOREIGN KEY (param_id) REFERENCES public.param(param_id);


--
-- Name: concr_meas ftk_meas; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_meas
    ADD CONSTRAINT ftk_meas FOREIGN KEY (meas_id) REFERENCES public.measurement(meas_id);


--
-- Name: concr_model ftk_model; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.concr_model
    ADD CONSTRAINT ftk_model FOREIGN KEY (mod_id) REFERENCES public.model(mod_id);


--
-- Name: model_condition_interpol ftkmod_id; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.model_condition_interpol
    ADD CONSTRAINT ftkmod_id FOREIGN KEY (mod_id) REFERENCES public.model(mod_id);


--
-- Name: model_condition_interpol ftkparam_id; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.model_condition_interpol
    ADD CONSTRAINT ftkparam_id FOREIGN KEY (param_id) REFERENCES public.param(param_id);


--
-- Name: measinput measinput_meas_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.measinput
    ADD CONSTRAINT measinput_meas_id_fkey FOREIGN KEY (meas_id) REFERENCES public.measurement(meas_id) ON UPDATE CASCADE;


--
-- Name: measinput measinput_param_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.measinput
    ADD CONSTRAINT measinput_param_id_fkey FOREIGN KEY (param_id) REFERENCES public.param(param_id) ON UPDATE CASCADE;


--
-- Name: measoutput measoutput_meas_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.measoutput
    ADD CONSTRAINT measoutput_meas_id_fkey FOREIGN KEY (meas_id) REFERENCES public.measurement(meas_id) ON UPDATE CASCADE;


--
-- Name: measoutput measoutput_param_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.measoutput
    ADD CONSTRAINT measoutput_param_id_fkey FOREIGN KEY (param_id) REFERENCES public.param(param_id) ON UPDATE CASCADE;


--
-- Name: model_condition_compare mod_id_1; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.model_condition_compare
    ADD CONSTRAINT mod_id_1 FOREIGN KEY (mod_id) REFERENCES public.model(mod_id);


--
-- Name: modelinput modelinput_mod_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modelinput
    ADD CONSTRAINT modelinput_mod_id_fkey FOREIGN KEY (mod_id) REFERENCES public.model(mod_id) ON UPDATE CASCADE;


--
-- Name: modelinput modelinput_param_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modelinput
    ADD CONSTRAINT modelinput_param_id_fkey FOREIGN KEY (param_id) REFERENCES public.param(param_id) ON UPDATE CASCADE;


--
-- Name: modeloutput modeloutput_mod_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modeloutput
    ADD CONSTRAINT modeloutput_mod_id_fkey FOREIGN KEY (mod_id) REFERENCES public.model(mod_id) ON UPDATE CASCADE;


--
-- Name: modeloutput modeloutput_param_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.modeloutput
    ADD CONSTRAINT modeloutput_param_id_fkey FOREIGN KEY (param_id) REFERENCES public.param(param_id) ON UPDATE CASCADE;


--
-- Name: param_semantic param_id; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.param_semantic
    ADD CONSTRAINT param_id FOREIGN KEY (param_id) REFERENCES public.param(param_id) ON UPDATE CASCADE;


--
-- Name: curve_params param_id; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.curve_params
    ADD CONSTRAINT param_id FOREIGN KEY (param_id) REFERENCES public.param(param_id);


--
-- Name: model_condition_compare param_id_1; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.model_condition_compare
    ADD CONSTRAINT param_id_1 FOREIGN KEY (param_id) REFERENCES public.param(param_id);


--
-- Name: sample sample_mat_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.sample
    ADD CONSTRAINT sample_mat_id_fkey FOREIGN KEY (mat_id) REFERENCES public.material(mat_id) ON UPDATE CASCADE;


--
-- Name: model typefk; Type: FK CONSTRAINT; Schema: public; Owner: mena
--

ALTER TABLE ONLY public.model
    ADD CONSTRAINT typefk FOREIGN KEY (mod_type) REFERENCES public.modeltype(modeltype) ON UPDATE CASCADE;


--
-- PostgreSQL database dump complete
--

