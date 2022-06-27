from openbabel import openbabel

def molecule_converter(molecule_input, input_format, output_format):
    """Converts chemical files given a file

    Parameters
    ----------
    molecule_input : str
        The file location of the chemical file
    input_format : str
        Format of the input file
    output_format : str
        Format of the output file

    Returns
    -------
    str
        a string with the content of the file, as well as the file itself
        
    -------
    Warning: openbabel somehow cannot convert gjf into xyz files
    """
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(input_format, output_format)
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, molecule_input)
    
    
    molecule_output = 'files' + '/' + molecule_input.split('/')[1].split('.')[0] + '.' +output_format
    obConversion.WriteFile(mol, molecule_output)
    mol_file = open(molecule_output)
    
    return mol_file.read()