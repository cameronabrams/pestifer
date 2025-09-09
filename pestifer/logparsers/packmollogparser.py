# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Packmol log parsing utility
"""

import logging
import os
import re
import yaml

import matplotlib.pyplot as plt
import pandas as pd

from .logparser import LogParser, get_single

logger = logging.getLogger(__name__)

class PackmolLogParser(LogParser):
    """
    A class for parsing Packmol log files. This class is a subclass of :class:`LogParser <pestifer.logparsers.logparser.LogParser>` and provides methods for reading, updating, and dumping Packmol log data.

    Parameters
    ----------
    basename : str
        The base name for the log parser. This is used to name the output log file.
    """
    banner_separator = '#' * 80  # banners can occur within sections
    """
    The separator used to identify banners in the Packmol log file.
    """
    section_separator = '-' * 80  # file is divided into sections
    """
    The separator used to identify sections in the Packmol log file.
    """
    def __init__(self, basename='packmol-logparser'):
        super().__init__()
        self.processed_separator_idx = []
        self.processed_sections = []
        self.processed_banners = []
        self.metadata = {}
        self.gencan = {}
        self.molmoves = []
        self.basename = basename
        self.progress = 0.0

    def update(self, bytes):
        """
        Update the Packmol log parser with new bytes of data. This method appends the new bytes to the byte collector and processes the sections in the log file.
        """
        super().update(bytes)
        separator_idx = [0] + [m.start() for m in re.finditer(self.section_separator, self.byte_collector)]
        for i, j in zip(separator_idx[:-1], separator_idx[1:]):
            if i not in self.processed_separator_idx:
                self.processed_separator_idx.append(i)
                self.process_section(self.byte_collector[i:j])

    def measure_progress(self):
        """
        Measure the progress of the Packmol log parsing. This method calculates the progress based on the number of GENCAN loops completed and the maximum number of GENCAN loops.
        It updates the `progress` attribute with the calculated progress value.
        """
        self.progress = self.metadata['current_total_gencan_loops'] / self.metadata['max_total_gencan_loops'] if 'current_total_gencan_loops' in self.metadata and 'max_total_gencan_loops' in self.metadata else 0.0
        return super().measure_progress()

    def process_banner(self, bytes: str):
        """
        Process a banner in the Packmol log file. This method identifies the type of banner and updates the metadata accordingly.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the banner to be processed. This can be a string or bytes object containing the banner data.
        """
        if 'Version' in bytes:
            vidx = bytes.index('Version') + len('Version') + 1
            vstr = bytes[vidx:].split()[0]
            self.metadata['version'] = vstr
            self.phase = 'initialization'
        elif 'Building initial approximation' in bytes:
            self.metadata['initial_approximation'] = {}
            self.phase = 'initial_approximation'
        elif 'Objective function at initial point' in bytes:
            idx = bytes.index('Objective function at initial point') + len('Objective function at initial point') + 1
            ofstr = bytes[idx:].split()[0]
            self.metadata['initial_objective_function'] = float(ofstr)
            self.phase = 'packing_molecules'
        elif 'Packing all molecules together' in bytes:
            self.phase = 'packing_all_molecules_together'
        elif 'Packing molecules of type' in bytes:
            idx = bytes.index('Packing molecules of type') + len('Packing molecules of type') + 1
            mtype = bytes[idx:].split()[0]
            self.phase = f'packing_molecules_of_type_{mtype}'
            if not 'molecule_types_packed' in self.metadata:
                self.metadata['molecule_types_packed'] = []
            self.metadata['molecule_types_packed'].append(mtype)
        self.processed_banners.append(bytes)

    def process_header(self, bytes: str):
        """
        Process the header section of the Packmol log file. This method extracts metadata from the header and updates the `metadata` attribute with relevant information.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the header section to be processed. This can be a string or bytes object
            containing the header data."""
        self.metadata['url'] = get_single('Userguide at:', bytes)
        flag = 'Types of coordinate files specified:'
        self.metadata['coordinate_file_types'] = bytes[bytes.index(flag) + len(flag):].split()[0]
        flag = 'Periodic boundary condition activated:'
        self.metadata['pbc_boxsize'] = [float(bytes[bytes.index(flag) + len(flag):].split()[x].strip()) for x in range(0, 3)]
        flag = 'PBC Reference box:'
        self.metadata['pbc_reference_box'] = [float(bytes[bytes.index(flag) + len(flag):].split()[x].strip()) for x in range(0, 6)]
        flag = 'Seed for random number generator:'
        self.metadata['seed'] = int(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
        flag = 'Output file:'
        self.metadata['output_file'] = bytes[bytes.index(flag) + len(flag):].split()[0].strip()
        cf_idx = [m.start() for m in re.finditer('Reading coordinate file:', bytes)]
        self.metadata['coordinate_files'] = []
        for i in cf_idx:
            idx = i + len('Reading coordinate file:')
            self.metadata['coordinate_files'].append(bytes[idx:].split()[0].strip())
        flag = 'Number of independent structures:'
        self.metadata['number_of_independent_structures'] = int(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
        str_idx = [m.start() for m in re.finditer('Structure', bytes)]
        # logger.debug(f'Found {len(str_idx)} structures')
        self.metadata['structures'] = []
        for i in str_idx:
            idx = i + len('Structure')
            eol = bytes[idx:].index('\n') + idx
            substr = bytes[idx:eol]
            # logger.debug(f'Found structure: {substr} {idx} {eol}')
            tokens = substr.replace('(', '').replace(')', '').replace(':', '').split()
            idx = int(tokens[0])
            fl = tokens[1]
            na = int(tokens[2])
            self.metadata['structures'].append(dict(idx=idx, file=fl, natoms=na))
        flag = 'Maximum number of GENCAN loops for all molecule packing:'
        self.metadata['max_gencan_loops_all'] = int(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
        flag = 'Maximum number of GENCAN loops for type:'
        mgl_idx = [m.start() for m in re.finditer(flag, bytes)]
        for i in mgl_idx:
            idx = i + len(flag)
            eol = bytes[idx:].index('\n') + idx
            substr = bytes[idx:eol]
            tokens = substr.replace(':', '').split()
            mtype = int(tokens[0])
            mloops = int(tokens[1])
            match = next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
            if match:
                match['max_gencan_loops'] = mloops
                match['current_gencan_loops'] = 0
            else:
                logger.error(f'process_header: No match for mtype {mtype} in structures {self.metadata["structures"]}')
        self.metadata['max_total_gencan_loops'] = self.metadata['max_gencan_loops_all'] + sum([d['max_gencan_loops'] for d in self.metadata['structures'] if 'max_gencan_loops' in d])
        self.metadata['current_total_gencan_loops'] = 0
        flag = 'Distance tolerance:'
        self.metadata['distance_tolerance'] = float(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
        flag = 'Number of molecules of type'
        nm_idx = [m.start() for m in re.finditer(flag, bytes)]
        for i in nm_idx:
            idx = i + len(flag)
            eol = bytes[idx:].index('\n') + idx
            substr = bytes[idx:eol]
            tokens = substr.replace(':', '').split()
            mtype = int(tokens[0])
            n = int(tokens[1])
            match = next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
            if match:
                match['number_of_molecules'] = n
        flag = 'Total number of atoms:'
        self.metadata['total_number_atoms'] = int(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
        flag = 'Total number of molecules:'
        self.metadata['total_number_molecules'] = int(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
        flag = 'Number of fixed molecules:'
        self.metadata['number_of_fixed_molecules'] = int(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
        flag = 'Number of free molecules:'
        self.metadata['number_of_free_molecules'] = int(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
        flag = 'Number of variables:'
        self.metadata['number_of_variables'] = int(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
        flag = 'Total number of fixed atoms:'
        self.metadata['total_number_fixed_atoms'] = int(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
        flag = 'Maximum internal distance of type'
        mid_idx = [m.start() for m in re.finditer(flag, bytes)]
        for i in mid_idx:
            idx = i + len(flag)
            eol = bytes[idx:].index('\n') + idx
            substr = bytes[idx:eol]
            tokens = substr.replace(':', '').split()
            mtype = int(tokens[0])
            mid = float(tokens[1])
            match = next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
            if match:
                match['maximum_internal_distance'] = mid

    def process_molecule_report(self, bytes: str):
        """
        Process a report about molecules in the Packmol log file. This method extracts information about the molecules, such as their type, function values before and after moving, and restraint violations.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the molecule report to be processed. This can be a string or bytes object
            containing the report data.
        """
        flag = 'Molecules of type'
        idx = bytes.index(flag) + len(flag)
        eol = bytes[idx:].index('\n') + idx
        substr = bytes[idx:eol]
        tokens = substr.replace(':', '').split()
        mtype = int(tokens[0])
        match = next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
        subflag = 'before_adjusting'
        if 'Adjusting random positions to fit the constraints.' in bytes:
            subflag = 'after_adjusting'
        # logger.debug(f'subflag: {subflag}')
        if match:
            match[subflag] = {}
            flag = 'Function value before moving molecules:'
            match[subflag]['function_value_before_moving_molecules'] = []
            fbm_idx = [m.start() for m in re.finditer(flag, bytes)]
            for i in fbm_idx:
                idx = i + len(flag)
                eol = bytes[idx:].index('\n') + idx
                substr = bytes[idx:eol]
                tokens = substr.replace(':', '').split()
                val = float(tokens[0])
                match[subflag]['function_value_before_moving_molecules'].append(val)
            flag = 'Function value after moving molecules:'
            match[subflag]['function_value_after_moving_molecules'] = []
            fam_idx = [m.start() for m in re.finditer(flag, bytes)]
            for i in fam_idx:
                idx = i + len(flag)
                eol = bytes[idx:].index('\n') + idx
                substr = bytes[idx:eol]
                tokens = substr.replace(':', '').split()
                val = float(tokens[0])
                match[subflag]['function_value_after_moving_molecules'].append(val)
            flag = 'Restraint-only function value:'
            match[subflag]['restraint_only_function_value'] = float(bytes[bytes.index(flag) + len(flag):].split()[0].strip())
            flag = 'Maximum violation of the restraints:'
            match[subflag]['maximum_violation_of_the_restraints'] = float(bytes[bytes.index(flag) + len(flag):].split()[0].strip())

    def process_gencan_report(self, bytes: str):
        """
        Process a report about the GENCAN (GENeric CANonical) packing algorithm in the Packmol log file. This method extracts information about the GENCAN loops, such as the iteration number, radius scaling, function values, and violations of target distances and constraints.
        It updates the ``gencan`` attribute with the extracted information.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the GENCAN report to be processed. This can be a string or bytes object
            containing the report data."""
        G = None
        ptok = self.phase.split('_')
        if ptok[1] == 'all':
            mtype = 'all'
        else:
            mtype = int(ptok[4])
        # logger.debug(f'process_gencan_report: {mtype}')
        if mtype not in self.gencan:
            self.gencan[mtype] = []
        G = self.gencan[mtype]
        if G is None:
            return
        iternum = int(get_single('Starting GENCAN loop:', bytes))
        radscal = float(get_single('Scaling radii by:', bytes))
        fval_last = float(get_single('Function value from last GENCAN loop: f =', bytes))
        fval_best = float(get_single('Best function value before: f =', bytes))
        max_viol_dist = float(get_single('Maximum violation of target distance:', bytes))
        max_viol_constr = float(get_single('Maximum violation of the constraints:', bytes))
        dsofar = dict(iteration=iternum, radius_scaling=radscal, function_value_last=fval_last, function_value_best=fval_best, max_viol_dist=max_viol_dist, max_viol_constr=max_viol_constr)
        if 'All-type function value:' in bytes:
            all_type_f = float(get_single('All-type function value:', bytes))
            dsofar['all_type_function_value'] = all_type_f
        G.append(dsofar)
        self.metadata['current_total_gencan_loops'] += 1
        match = next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
        if match:
            match['current_gencan_loops'] += 1

    def process_gencan_success(self, bytes: str):
        """
        Process a report about the successful completion of the GENCAN packing algorithm in the Packmol log file. This method extracts information about the molecule type, objective function value, maximum violation of target distance, and maximum constraint violation.
        It updates the metadata of the corresponding molecule type in the ``metadata['structures']`` attribute
        with the extracted information.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the GENCAN success report to be processed. This can be a string or bytes object
            containing the report data. 
        """
        mtype = int(get_single('Packing solved for molecules of type', bytes))
        obj_func_val = float(get_single('Objective function value:', bytes))
        max_viol_target_dist = float(get_single('Maximum violation of target distance:', bytes))
        max_viol_constr = float(get_single('Max. constraint violation:', bytes))
        match = next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
        if match:
            match['gencan_success'] = dict(objective_function_value=obj_func_val, max_viol_target_distance=max_viol_target_dist, max_viol_constr=max_viol_constr)
            unneeded_loops = match.get('max_gencan_loops', 0) - match.get('current_gencan_loops', 0)
            self.metadata['max_total_gencan_loops'] -= (unneeded_loops - 1)

    def process_moving_worst_molecules(self, bytes: str):
        """
        Process a report about moving the worst molecules in the Packmol log file. This method extracts information about the function values before and after moving molecules, the types of molecules moved, and their percentages.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the report about moving worst molecules to be processed. This can be a string or bytes object
            containing the report data.
        """
        top_idx = [m.start() for m in re.finditer('Moving', bytes)]
        cont_idx = top_idx[0] + len('Moving worst molecules ...')
        mbytes = bytes[cont_idx:]
        fvbefore = float(get_single('Function value before moving molecules:', mbytes))
        flag = 'Moving '
        moltypes = []
        pct = []
        fam_idx = [m.start() for m in re.finditer(flag, mbytes)]
        for i in fam_idx:
            idx = i + len(flag)
            eol = mbytes[idx:].index('\n') + idx
            substr = mbytes[idx:eol]
            tokens = substr.replace(':', '').split()
            val = int(tokens[0])
            cnt = int(tokens[-1])
            moltypes.append((val, cnt))
        flag = 'Type'
        fam_idx = [m.start() for m in re.finditer(flag, mbytes)]
        for i in fam_idx:
            idx = i + len(flag)
            eol = mbytes[idx:].index('\n') + idx
            substr = mbytes[idx:eol]
            tokens = substr.replace(':', '').split()
            val = int(tokens[0])
            cnt = float(tokens[-1].replace('%', ''))
            pct.append((val, cnt))
        fvafter = float(get_single('Function value after moving molecules:', mbytes))
        result = dict(function_value_before_moving_molecules=fvbefore,
                      function_value_after_moving_molecules=fvafter,
                      moltypes=moltypes,
                      pcts=pct)
        self.molmoves.append(result)

    def process_section(self, bytes: str):
        """
        Process a section of the Packmol log file. This method identifies the type of section and processes it accordingly.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the section to be processed. This can be a string or bytes object containing
            the section data.
        """
        banner_idx = [m.start() for m in re.finditer(self.banner_separator, bytes)]
        nbsep = len(banner_idx)
        banners = []
        if nbsep > 1:
            for i, j in zip(banner_idx[:-1:2], banner_idx[1::2]):
                banners.append(bytes[i:j])
            if nbsep % 2 == 1:
                banners.append(bytes[banner_idx[-2]:banner_idx[-1]])
        for b in banners:
            self.process_banner(b)

        if 'Reading input file...' in bytes:
            self.process_header(bytes)
        elif 'Molecules of type:' in bytes:
            self.process_molecule_report(bytes)
        elif 'Starting GENCAN loop:' in bytes:
            self.process_gencan_report(bytes)
        elif 'Packing solved for molecules of type' in bytes:
            self.process_gencan_success(bytes)
        elif 'Moving worst molecules' in bytes:
            self.process_moving_worst_molecules(bytes)

    def finalize(self):
        """
        Finalize the Packmol log parser by saving the metadata and GENCAN data to files. This method writes the metadata to a YAML file and the GENCAN data to CSV files.
        It also generates plots of the GENCAN function values over iterations for each molecule type and saves them as PNG files.
        
        Returns
        -------
        str
            The filename of the generated PNG file containing the GENCAN function value plots.  
        """
        with open(f'{self.basename}_packmol-results.yaml', 'w') as f:
            yaml.dump(self.metadata, f, default_flow_style=False)
        self.gencan_df = {}
        for k, v in self.gencan.items():
            self.gencan_df[k] = pd.DataFrame(v)
            self.gencan_df[k].to_csv(f'{self.basename}_{k}_packmol.csv', index=False)
        fig, ax = plt.subplots(1, len(self.gencan_df), figsize=(4 * len(self.gencan_df), 4))
        for i, (k, v) in enumerate(self.gencan_df.items()):
            if len(v) == 0:
                continue
            ax[i].plot(v['iteration'], v['function_value_last'], label='Last function value')
            ax[i].plot(v['iteration'], v['function_value_best'], label='Best function value')
            if k == 'all':
                title = 'All molecules'
            else:
                assert isinstance(k, int), f'Expected k to be int, got {type(k)}'
                pdb_name = list(filter(lambda x: x['idx'] == k, self.metadata['structures']))[0]['file']
                mname = os.path.splitext(os.path.basename(pdb_name))[0]
                title = f'Molecule {mname} ({k})'
            ax[i].set_title(title)
            ax[i].set_xlabel('Iteration')
            ax[i].set_ylabel('Function value')
            ax[i].legend()
            ax[i].set_yscale('log')
            ax[i].grid(True)
        plt.tight_layout()
        plt.savefig(f'{self.basename}_packmol.png')
        plt.close()
        return f'{self.basename}_packmol.png'
