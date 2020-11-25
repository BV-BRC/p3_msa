package com.ngc.brc.analysisservices.snp;

import java.io.File;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import com.ngc.brc.commons.constants.Constants;
import com.ngc.brc.commons.utils.io.IoUtils;

public class SNPRunner implements Runnable {

	public static final String DIRECTORY_FOR_SNP_FILES = Constants.DIRECTORY_FOR_ANALYSIS_FILES;
	public static final String FINAL_AFA_FILE_NAME = "output.afa";
	public static final String FINAL_OUTPUT_FILE_NAME = "foma.table";
	private static Log logger = LogFactory.getLog(SNPRunner.class);

	private String ticketNumber;
	private String options;
	private String sequence;
	private String parameters;

	public void initialize(String ticketNumber, String options, String parameters, String sequence) {
		this.ticketNumber = ticketNumber;
		this.options = options;
		this.parameters = parameters;
		this.sequence = sequence;
	}

	public void run() {

		File folder = new File(DIRECTORY_FOR_SNP_FILES + File.separator + ticketNumber);
		if (!folder.exists()) {
			folder.mkdir();
		}
		if (options.lastIndexOf("-x") > 0) {
			IoUtils.writeStringToFile(sequence, DIRECTORY_FOR_SNP_FILES + File.separator + ticketNumber + File.separator + FINAL_AFA_FILE_NAME);
		} else {
			IoUtils.writeStringToFile(sequence, DIRECTORY_FOR_SNP_FILES + File.separator + ticketNumber + File.separator + Constants.FASTA_INPUT_FILE_NAME);
		}
		IoUtils.writeStringToFile(parameters, DIRECTORY_FOR_SNP_FILES + File.separator + ticketNumber + File.separator + Constants.PARAM_FILE_NAME);
		String cmd = "web_flu_snp_analysis.pl -r " + DIRECTORY_FOR_SNP_FILES + File.separator + ticketNumber + " " + options;
		logger.info("SNP " + cmd);
		IoUtils.runCommand(cmd, null, null);
	}
}