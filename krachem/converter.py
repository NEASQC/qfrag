from openbabel import openbabel

def molecule_converter(molecule_input,input_format, output_format):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(input_format, output_format)
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, molecule_input)
    
    
    molecule_output = 'files' + '/' + molecule_input.split('/')[1].split('.')[0] + '.' +output_format
    obConversion.WriteFile(mol, molecule_output)
    mol_file = open(molecule_output)
    
    return mol_file.read()

def gjf2xyz(molecule_input):
    molecule_content = open(molecule_input,'r').readlines()
    molecule_content[0] = (molecule_content[-1].strip()) + '\n'
    del molecule_content[2:6]
    del molecule_content[(int(molecule_content[0])+2)::]
    
    molecule_output = 'xyz' + '/' + molecule_input.split('/')[1].split('.')[0]
    
    with open(molecule_output+'.xyz', "w") as file:
        [file.write(x) for x in molecule_content]
       
    mol_file = open(molecule_output+'.xyz')
    return mol_file.read()