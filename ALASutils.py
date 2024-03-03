from astropy.table import Table
import pandas as pd

from functools import reduce

def ReadCatalogue ( file, log = False ):
    """ SourceExtractor y SourceXtractor++ tiene salidas con columnas de diferentes dimensiones: 2 o 3. 
    Por eso se usa el condicional para separar las columnas dependiendo del t.shape. 
    
    Modo de uso: catalogue = ReadCatalogue ( file ) """
    tmp = pd.DataFrame()
    
    t = Table.read ( file, format = "fits" )
    print ("[*] Cantidad de datos: %i (%s)" % (len(t), file) )

    for j, c in enumerate(t.colnames): # For multidimensional columns
        if len(t[c].shape) > 1: 
            for col in [ c ]:
                if log: print ( c, len(t[c].shape), t[c].shape )
                for i in range(len(t[c].shape)):
                    if (len(t[c].shape) == 2):
                        t [col + str(i)] = t [col][:,i].flatten()
                    elif (len(t[c].shape) == 3) & (t[c].shape[1] == 1):
                        t [col + str(i)] = t [col][:,:,i].flatten()
                    elif (len(t[c].shape) == 3) & (t[c].shape[1] > 1): 
                        t [col + str(i)] = t [col][:,i].flatten()
                    else:
                        print ( '[!] Falta definir algun tipo mas de columna:', t[c].shape )
                        return
                        
    names = [ name for name in t.colnames if len(t[name].shape) <= 1 ] 
    tmp = t[names].to_pandas()
    return tmp

def DownloadCatalogue (url, readme):
    """ Descarga catalogos de CDS y los convierte a pandas. 
    Es necesario indicar el archivo README para que interprete el formato.
    
    Modo de uso: catalogue = DownloadCatalogue (url, readme) """
    from astropy.io import ascii
    return ascii.read ( url, readme = readme ).to_pandas()

def MergeCatalogues(frames, bands):
    """ Une dos o mas catalogos. Toma los primeros dos catalogos para identificar las columnas en comun. Luego une todos los
    catalogos, agregando las columnas en comun una sola vez.
    Esto fue pensado para usar con catalogos salidos de SourcExtractor.

    MODO DE USO: df = MergeCatalogues( [catR, catG, catZ], ['R', 'G', 'Z'] )    
    """
    if len(frames) != len(bands): raise ValueError('Dimension de frames y bandas es distinta.')
    cat = pd.DataFrame()
    
    if len(frames) >= 2:
        identical_columns = [ column for column in frames[0].columns if frames[0][column].equals(frames[1][column]) ]
        new_frames = []
        for i, frame in enumerate(frames):
            new_frame = frame.drop(columns = identical_columns)
            new_frame.columns = [ c + '_' + str(bands[i]) for c in new_frame.columns ]
            new_frames.append(new_frame)

        cat = frames[0][identical_columns]
        return reduce(lambda left, right: pd.merge(left,right,left_index = True, right_index = True), [cat] + new_frames)
        
    else:
        raise ValueError('Necesito 2 o mas catalogos.')
        
def CrossMatch(cat1, cat2, coor1 = ['RA', 'DEC'], coor2 = ['RA', 'DEC'], sep_min = 1):
    """ Cross-match entre dos catalogos. La funcion no limita el resultado de cada coincidencia.

    Modo de uso: CM = CrossMatch(cat1, cat2, coor1 = ['RA', 'DEC'], coor2 = ['RA', 'DEC'], sep_min = 3 * u.arcsec) """
    c2 = SkyCoord(cat2[coor2[0]], cat2[coor2[1]], unit = (u.deg, u.deg))
    c1 = SkyCoord(cat1[coor1[0]], cat1[coor1[1]], unit = (u.deg, u.deg))

    idx_c2, idx_c1, d2d, d3d = c1.search_around_sky ( c2, sep_min, storekdtree = True )

    cat2_match = cat2.iloc[ idx_c2 ]
    cat2_match.reset_index ( inplace = True )
    cat2_match = cat2_match.add_suffix('_c2')
    
    cat1_match = cat1.iloc[ idx_c1 ]
    cat1_match.reset_index ( inplace = True )
    cat1_match = cat1_match.add_suffix('_c1')

    Xmatch = pd.concat ( [cat1_match, cat2_match], axis = 1 )
    return Xmatch

# Ver si esto se puede evitar, usando la columna d2d que sale del search_around_sky
# SkyCoord ( Xmatch[['RA_c1','DEC_c1']].values, unit = (u.deg, u.deg)).separation( SkyCoord( Xmatch[['RA_c2','DEC_c2']].values, unit = (u.deg, u.deg)) ).arcsec
# Xmatch.sort_values ( 'separation' )
# Xmatch.drop_duplicates ( subset = 'index_c2', keep = 'first' ).drop_duplicates ( subset = 'index_c1', keep = 'first' )
