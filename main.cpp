#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>


struct process_Info {
	char processID;
	int arrivalTime;
	int numBursts;
	std::vector<int> cpuBurstTimes;
	std::vector<int> ioBurstTimes;
	std::vector<int> turnaroundTimes;
	std::vector<int> waitTimes;
	int burstsCompleted = 0;
	int burstStartTime = -1;
	int burstTimeRemaining = 0;
	int burstEndTime = -1;

	//returns true if this should come before p2 in an SJF algorithm
	bool compareSJF(process_Info& p2) { //TODO: make this fit the pdf specifications better
		if (cpuBurstTimes[burstsCompleted] < p2.cpuBurstTimes[p2.burstsCompleted]) {
			return true;
		} else if (cpuBurstTimes[burstsCompleted] < p2.cpuBurstTimes[p2.burstsCompleted]) {
			if (processID < p2.processID) {
				return true;
			}
		}
		return false;
	}
};

double next_exp(int max, double lambda);
void FCFS(std::vector<process_Info> processes, int csTime);
void SJF(std::vector<process_Info> processes, int csTime);
void SRT(std::vector<process_Info> processes, int csTime);
void RR(std::vector<process_Info> processes, int csTime, int timeSlice);

int main(int argc, char ** argv){
	//Error Checking Arguments
	if(argc != 8){
		fprintf(stderr, "Not enough command-line arguments.\n");
		return EXIT_FAILURE;
	}
	int numProc = atoi(argv[1]);
	if(numProc < 1 || numProc > 26){
		fprintf(stderr, "Invalid number of processes.\n");
		return EXIT_FAILURE;
	}
	int seed = atoi(argv[2]);
	if (seed == 0) {
		fprintf(stderr, "Invalid seed.\n");
		return EXIT_FAILURE;
	}
	double lambda = atof(argv[3]);
	if (lambda == 0.0f) {
		fprintf(stderr, "Invalid lambda.\n");
		return EXIT_FAILURE;
	}
	int upBound = atoi(argv[4]);
	if (upBound == 0) {
		fprintf(stderr, "Invalid upper bound.\n");
		return EXIT_FAILURE;
	}
	int csTime = atoi(argv[5]);
	if (csTime == 0) {
		fprintf(stderr, "Invalid context switch time.\n");
		return EXIT_FAILURE;
	}
	double alpha = atof(argv[6]);
	if (alpha == 0.0f) {
		fprintf(stderr, "Invalid alpha.\n");
		return EXIT_FAILURE;
	}
	int timeSlice = atoi(argv[7]);
	if (timeSlice == 0) {
		fprintf(stderr, "Invalid time slice.\n");
		return EXIT_FAILURE;
	}
	//Setting up information for each process
	srand48(seed);
	std::vector<process_Info> processes;
	for (int i = 0; i < numProc; i++) {
		process_Info process;
		process.processID = 65 + i;
		process.arrivalTime = floor(next_exp(upBound, lambda));
		process.numBursts = ceil(drand48() * 100);
		process.waitTimes.resize(process.numBursts);
		process.turnaroundTimes.resize(process.numBursts);
		for (int j = 0; j < process.numBursts; j++) {
			process.cpuBurstTimes.push_back(ceil(next_exp(upBound, lambda)));
			if (j != process.numBursts - 1) {
				process.ioBurstTimes.push_back(ceil(next_exp(upBound, lambda)) * 10);
			}
		}
		processes.push_back(process);
		std::string processInfo = "Process ";
		processInfo += std::string(1, process.processID);
		processInfo += " (arrival time ";
		processInfo += std::to_string(process.arrivalTime);
		processInfo += " ms) ";
		processInfo += std::to_string(process.numBursts);
		processInfo += " CPU bursts (tau ";
		processInfo += std::to_string(int(ceil(1 / lambda)));
		processInfo += "ms)";
		std::cout << processInfo << std::endl;
	}
	std::cout << std::endl;
	FCFS(processes, csTime);
	SJF(processes, csTime);
	RR(processes, csTime, timeSlice);
}

//Calculating a random number based on exponential distribution
double next_exp(int max, double lambda) {
	double x = 0;
	while (true)
	{
		double r = drand48();   /* uniform dist [0.00,1.00) -- also see random() */

		x = -log(r) / lambda;  /* log() is natural log */

		/* skip values that are far down the "long tail" of the distribution */
		if (x > max) { continue; }

		break;
	}
	return x;
}

//Returns the ready queue in the format [Q [List of processes in the ready queue]]
//Or [Q empty] if the ready queue is empty
std::string ready_Queue_Format(std::list<process_Info*> readyQueue) {
	std::string queue = "[Q ";
	if (readyQueue.size() == 0) {
		queue += "empty";
	}
	else {
		for (process_Info * p : readyQueue) {
			queue += p->processID;
		}
	}
	queue += "]";
	return queue;
}

void FCFS(std::vector<process_Info> processes, int csTime) {
	std::cout << "time 0ms: Simulator started for FCFS [Q empty]" << std::endl;
	int time = 0;
	std::list<process_Info*> readyQueue;
	process_Info * runningProcess;
	bool running = false;
	int cpuEmptyTime = 0;
	float cpuUtilization = 0;
	int numContextSwitches = 0;
	bool cpuUtilized = false;
	while (true) {
		//Context switching the process out of the CPU is done
		if (cpuEmptyTime == time) {
			running = false;
			runningProcess = nullptr;
			bool done = true;
			for (int i = 0; i < int(processes.size()); i++) {
				if (processes[i].burstEndTime != 0) {
					done = false;
				}
			}
			if (done) {
				std::cout << "time " << time << "ms: Simulator ended for FCFS " << ready_Queue_Format(readyQueue) << std::endl;
				break;
			}
		}
		//CPU is being utilized(not including context switches)
		if (cpuUtilized)
			cpuUtilization++;
		//Running process switched in ( context switch done )
		if (running && runningProcess->burstStartTime == time) {
			cpuUtilized = true;
			runningProcess->turnaroundTimes[runningProcess->burstsCompleted] += csTime / 2;
			std::cout << "time " << time << "ms: Process " << runningProcess->processID << " started using the CPU for " 
					  << runningProcess->cpuBurstTimes[runningProcess->burstsCompleted] << "ms burst " << ready_Queue_Format(readyQueue) << std::endl;
		}
		//Process in the CPU done running, still need to account for context switch times
		if (running && runningProcess->burstEndTime == time) {
			cpuUtilized = false;
			cpuEmptyTime = time + csTime / 2;
			numContextSwitches++;
			runningProcess->turnaroundTimes[runningProcess->burstsCompleted] += runningProcess->burstEndTime - runningProcess->burstStartTime + csTime / 2;
			runningProcess->burstsCompleted++;
			std::string burst = "";
			//If this is the last burst for this process
			if (runningProcess->burstsCompleted == int(runningProcess->cpuBurstTimes.size())) {
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " terminated "
						  << ready_Queue_Format(readyQueue) << std::endl;
				runningProcess->burstEndTime = 0;
			}
			else {
				if (runningProcess->burstsCompleted == int(runningProcess->cpuBurstTimes.size()) - 1)
					burst = "burst";
				else
					burst = "bursts";
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " completed a CPU burst; "
					<< runningProcess->cpuBurstTimes.size() - runningProcess->burstsCompleted << " " << burst << " to go " << ready_Queue_Format(readyQueue) << std::endl;
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " switching out of CPU; will block on I/O until time "
					<< time + runningProcess->ioBurstTimes[runningProcess->burstsCompleted - 1] + csTime / 2 << "ms " << ready_Queue_Format(readyQueue) << std::endl;
				runningProcess->burstEndTime = time + runningProcess->ioBurstTimes[runningProcess->burstsCompleted - 1] + csTime / 2;
			}
		}
		//Checking if each process has completed its I/O burst time
		for (int i = 0; i < int(processes.size()); i++) {
			if (processes[i].burstEndTime == time) {
				readyQueue.push_back(&processes[i]);
				std::cout << "time " << time << "ms: Process " << processes[i].processID << " completed I/O; added to ready queue " << ready_Queue_Format(readyQueue) << std::endl;
			}
		}
		//Checking if each process has arrived
		for (int i = 0; i < int(processes.size()); i++) {
			if (processes[i].arrivalTime == time) {
				readyQueue.push_back(&processes[i]);
				std::cout << "time " << time << "ms: Process " << processes[i].processID << " arrived; added to ready queue " << ready_Queue_Format(readyQueue) << std::endl;
			}
		}
		//If there's no process running, wait for the context switch time and set the start and end burst times
		if (!running && readyQueue.size() > 0 && (*readyQueue.begin())->burstStartTime < time) {
			runningProcess = (*readyQueue.begin());
			running = true;
			readyQueue.pop_front();
			runningProcess->burstStartTime = time + csTime / 2;
			runningProcess->burstEndTime = runningProcess->burstStartTime + runningProcess->cpuBurstTimes[runningProcess->burstsCompleted];
		}
		//Adding to each process's wait time and turnaround time when in the ready queue
		for (std::list<process_Info*>::iterator itr = readyQueue.begin(); itr != readyQueue.end(); itr++) {
			(*itr)->waitTimes[(*itr)->burstsCompleted]++;
			(*itr)->turnaroundTimes[(*itr)->burstsCompleted]++;
		}
		time++;
	}
	//Calculate all the statistics
	float avgWaitTime = 0;
	float avgturnAround = 0;
	float avgBurstTime = 0;
	int totalBursts = 0;
	for (int i = 0; i < int(processes.size()); i++) {
		for (int j = 0; j < processes[i].numBursts; j++) {
			avgturnAround += processes[i].turnaroundTimes[j];
			avgWaitTime += processes[i].waitTimes[j];
			avgBurstTime += processes[i].cpuBurstTimes[j];
			totalBursts++;
		}
	}
	avgturnAround /= totalBursts;
	avgWaitTime /= totalBursts;
	avgBurstTime /= totalBursts;
	std::ofstream simout;
	simout.open("simout.txt");
	simout << "Algorithm FCFS\n";
	simout << "-- average CPU burst time: " << std::fixed << std::setprecision(3) << avgBurstTime << " ms\n";
	simout << "-- average wait time: " << std::fixed << std::setprecision(3) << avgWaitTime << " ms\n";
	simout << "-- average turnaround time: " << std::fixed << std::setprecision(3) << avgturnAround << " ms\n";
	simout << "-- total number of context switches: " << numContextSwitches << std::endl;
	simout << "-- total number of preemptions: 0\n";
	simout << "-- CPU utilization: " << std::fixed << std::setprecision(3) << cpuUtilization / time * 100 << "%\n";
	simout.close();
}

void RR(std::vector<process_Info> processes, int csTime, int timeSlice) {
	std::cout << "time 0ms: Simulator started for RR with time slice " << timeSlice << "ms [Q empty]" << std::endl;
	int time = 0;
	std::list<process_Info*> readyQueue;
	process_Info* runningProcess;
	bool running = false;
	int cpuEmptyTime = -1;
	float cpuUtilization = 0;
	int numContextSwitches = 0;
	int timeSliceTime = -1;
	bool cpuUtilized = false;
	int numPreemptions = 0;
	while (true) {
		//Context switching the process out of the CPU is done
		if (cpuEmptyTime == time) {
			running = false;
			if (runningProcess->burstTimeRemaining > 0) {
				readyQueue.push_back(runningProcess);
			}
			runningProcess = nullptr;
			bool done = true;
			for (int i = 0; i < int(processes.size()); i++) {
				if (processes[i].burstEndTime != 0) {
					done = false;
				}
			}
			if (done) {
				std::cout << "time " << time << "ms: Simulator ended for RR " << ready_Queue_Format(readyQueue) << std::endl;
				break;
			}
		}
		//CPU is being utilized(not including context switches)
		if (cpuUtilized)
			cpuUtilization++;
		//Running process switched in ( context switch done )
		if (running && runningProcess->burstStartTime == time) {
			cpuUtilized = true;
			timeSliceTime = time + timeSlice;
			runningProcess->turnaroundTimes[runningProcess->burstsCompleted] += csTime / 2;
			//Check if the burst is a new one or old one that hasn't finished
			if (runningProcess->burstTimeRemaining > 0) {
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " started using the CPU for remaining "
					<< runningProcess->burstTimeRemaining << "ms of " << runningProcess->cpuBurstTimes[runningProcess->burstsCompleted] << "ms burst " 
					<< ready_Queue_Format(readyQueue) << std::endl;
			}
			else {
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " started using the CPU for "
					<< runningProcess->cpuBurstTimes[runningProcess->burstsCompleted] << "ms burst " << ready_Queue_Format(readyQueue) << std::endl;
			}
		}
		//Process in the CPU is done running, still need to account for context switch times
		if (running && runningProcess->burstEndTime == time) {
			cpuUtilized = false;
			cpuEmptyTime = time + csTime / 2;
			numContextSwitches++;
			runningProcess->burstTimeRemaining = 0;
			runningProcess->turnaroundTimes[runningProcess->burstsCompleted] += runningProcess->burstEndTime - runningProcess->burstStartTime + csTime / 2;
			runningProcess->burstsCompleted++;
			std::string burst = "";
			//If this is the last burst for this process
			if (runningProcess->burstsCompleted == int(runningProcess->cpuBurstTimes.size())) {
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " terminated "
					<< ready_Queue_Format(readyQueue) << std::endl;
				runningProcess->burstEndTime = 0;
			}
			else {
				if (runningProcess->burstsCompleted == int(runningProcess->cpuBurstTimes.size()) - 1)
					burst = "burst";
				else
					burst = "bursts";
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " completed a CPU burst; "
					<< runningProcess->cpuBurstTimes.size() - runningProcess->burstsCompleted << " " << burst << " to go " << ready_Queue_Format(readyQueue) << std::endl;
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " switching out of CPU; will block on I/O until time "
					<< time + runningProcess->ioBurstTimes[runningProcess->burstsCompleted - 1] + csTime / 2 << "ms " << ready_Queue_Format(readyQueue) << std::endl;
				runningProcess->burstEndTime = time + runningProcess->ioBurstTimes[runningProcess->burstsCompleted - 1] + csTime / 2;
			}
		}
		//If the timeslice is over
		if (cpuUtilized && timeSliceTime == time) {
			//Don't preempt if the ready queue is empty
			if (readyQueue.size() == 0) {
				std::cout << "time " << time << "ms: Time slice expired; no preemption because ready queue is empty " << ready_Queue_Format(readyQueue) << std::endl;
				timeSliceTime += timeSlice;
			}
			else{
				cpuUtilized = false;
				cpuEmptyTime = time + csTime / 2;
				numContextSwitches++;
				runningProcess->turnaroundTimes[runningProcess->burstsCompleted] += time - runningProcess->burstStartTime + csTime / 2;
				runningProcess->burstTimeRemaining = runningProcess->burstEndTime - time;
				runningProcess->burstEndTime = time - 1;
				std::cout << "time " << time << "ms: Time slice expired; process " << runningProcess->processID << " preempted with " 
						  << runningProcess->burstTimeRemaining << "ms to go " << ready_Queue_Format(readyQueue) << std::endl;
				numPreemptions++;
			}
		}
		//Checking if each process has completed its I/O burst time
		for (int i = 0; i < int(processes.size()); i++) {
			if (processes[i].burstEndTime == time) {
				readyQueue.push_back(&processes[i]);
				std::cout << "time " << time << "ms: Process " << processes[i].processID << " completed I/O; added to ready queue " << ready_Queue_Format(readyQueue) << std::endl;
			}
		}
		//Checking if each process has arrived
		for (int i = 0; i < int(processes.size()); i++) {
			if (processes[i].arrivalTime == time) {
				readyQueue.push_back(&processes[i]);
				std::cout << "time " << time << "ms: Process " << processes[i].processID << " arrived; added to ready queue " << ready_Queue_Format(readyQueue) << std::endl;
			}
		}
		//If there's no process running, wait for the context switch time and set the start and end burst times
		if (!running && readyQueue.size() > 0 && (*readyQueue.begin())->burstStartTime < time) {
			(*readyQueue.begin())->burstStartTime = time + csTime / 2;
			runningProcess = (*readyQueue.begin());
			readyQueue.pop_front();
			running = true;
			if (runningProcess->burstTimeRemaining > 0) {
				runningProcess->burstEndTime = runningProcess->burstStartTime + runningProcess->burstTimeRemaining;
			}
			else {
				runningProcess->burstEndTime = runningProcess->burstStartTime + runningProcess->cpuBurstTimes[runningProcess->burstsCompleted];
			}
		}
		//Adding to each process's wait time and turnaround time when in the ready queue
		for (std::list<process_Info*>::iterator itr = readyQueue.begin(); itr != readyQueue.end(); itr++) {
			(*itr)->waitTimes[(*itr)->burstsCompleted]++;
			(*itr)->turnaroundTimes[(*itr)->burstsCompleted]++;
		}
		time++;
	}
	//Calculate all the statistics
	float avgWaitTime = 0;
	float avgturnAround = 0;
	float avgBurstTime = 0;
	int totalBursts = 0;
	for (int i = 0; i < int(processes.size()); i++) {
		for (int j = 0; j < processes[i].numBursts; j++) {
			avgturnAround += processes[i].turnaroundTimes[j];
			avgWaitTime += processes[i].waitTimes[j];
			avgBurstTime += processes[i].cpuBurstTimes[j];
			totalBursts++;
		}
	}
	avgturnAround /= totalBursts;
	avgWaitTime /= totalBursts;
	avgBurstTime /= totalBursts;
	std::ofstream simout;
	simout.open("simout.txt", std::ios_base::app);
	simout << "Algorithm RR\n";
	simout << "-- average CPU burst time: " << std::fixed << std::setprecision(3) << avgBurstTime << " ms\n";
	simout << "-- average wait time: " << std::fixed << std::setprecision(3) << avgWaitTime << " ms\n";
	simout << "-- average turnaround time: " << std::fixed << std::setprecision(3) << avgturnAround << " ms\n";
	simout << "-- total number of context switches: " << numContextSwitches << std::endl;
	simout << "-- total number of preemptions: " << numPreemptions << std::endl;
	simout << "-- CPU utilization: " << std::fixed << std::setprecision(3) << cpuUtilization / time * 100 << "%\n";
	simout.close();
}

void SJF(std::vector<process_Info> processes, int csTime) {
	std::cout << "time 0ms: Simulator started for SJF [Q empty]" << std::endl;
	int time = 0;
	std::list<process_Info*> readyQueue;
	process_Info* runningProcess;
	bool running = false;
	int cpuEmptyTime = -1;
	float cpuUtilization = 0;
	int numContextSwitches = 0;
	int timeSliceTime = -1;
	bool cpuUtilized = false;
	int numPreemptions = 0;



while (true) {
		//Context switching the process out of the CPU is done
		if (cpuEmptyTime == time) {
			running = false;
			runningProcess = nullptr;
			bool done = true;
			for (int i = 0; i < int(processes.size()); i++) {
				if (processes[i].burstEndTime != 0) {
					done = false;
				}
			}
			if (done) {
				std::cout << "time " << time << "ms: Simulator ended for FCFS " << ready_Queue_Format(readyQueue) << std::endl;
				break;
			}
		}
		//CPU is being utilized(not including context switches)
		if (cpuUtilized)
			cpuUtilization++;
		//Running process switched in ( context switch done )
		if (running && runningProcess->burstStartTime == time) {
			cpuUtilized = true;
			runningProcess->turnaroundTimes[runningProcess->burstsCompleted] += csTime / 2;
			std::cout << "time " << time << "ms: Process " << runningProcess->processID << " started using the CPU for " 
					  << runningProcess->cpuBurstTimes[runningProcess->burstsCompleted] << "ms burst " << ready_Queue_Format(readyQueue) << std::endl;
		}
		//Process in the CPU done running, still need to account for context switch times
		if (running && runningProcess->burstEndTime == time) {
			cpuUtilized = false;
			cpuEmptyTime = time + csTime / 2;
			numContextSwitches++;
			runningProcess->turnaroundTimes[runningProcess->burstsCompleted] += runningProcess->burstEndTime - runningProcess->burstStartTime + csTime / 2;
			runningProcess->burstsCompleted++;
			std::string burst = "";
			//If this is the last burst for this process
			if (runningProcess->burstsCompleted == int(runningProcess->cpuBurstTimes.size())) {
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " terminated "
						  << ready_Queue_Format(readyQueue) << std::endl;
				runningProcess->burstEndTime = 0;
			}
			else {
				if (runningProcess->burstsCompleted == int(runningProcess->cpuBurstTimes.size()) - 1)
					burst = "burst";
				else
					burst = "bursts";
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " completed a CPU burst; "
					<< runningProcess->cpuBurstTimes.size() - runningProcess->burstsCompleted << " " << burst << " to go " << ready_Queue_Format(readyQueue) << std::endl;
				std::cout << "time " << time << "ms: Process " << runningProcess->processID << " switching out of CPU; will block on I/O until time "
					<< time + runningProcess->ioBurstTimes[runningProcess->burstsCompleted - 1] + csTime / 2 << "ms " << ready_Queue_Format(readyQueue) << std::endl;
				runningProcess->burstEndTime = time + runningProcess->ioBurstTimes[runningProcess->burstsCompleted - 1] + csTime / 2;
			}
		}
		//Checking if each process has completed its I/O burst time
		for (int i = 0; i < int(processes.size()); i++) {
			if (processes[i].burstEndTime == time) {
				//find spot where process has shorter runtime
				for (std::list<process_Info*>::iterator itr = readyQueue.begin(); itr != readyQueue.end(); itr++) { 
					if (processes[i].compareSJF(**itr)) {
						readyQueue.insert(itr, &processes[i]);
						break;
					}
				}
				std::cout << "time " << time << "ms: Process " << processes[i].processID << " completed I/O; added to ready queue " << ready_Queue_Format(readyQueue) << std::endl;
			}
		}
		//Checking if each process has arrived
		for (int i = 0; i < int(processes.size()); i++) {
			if (processes[i].arrivalTime == time) {
				//find spot where process has shorter runtime
				for (std::list<process_Info*>::iterator itr = readyQueue.begin(); itr != readyQueue.end(); itr++) { 
					if (processes[i].compareSJF(**itr)) {
						readyQueue.insert(itr, &processes[i]);
						break;
					}
				}
				std::cout << "time " << time << "ms: Process " << processes[i].processID << " arrived; added to ready queue " << ready_Queue_Format(readyQueue) << std::endl;
			}
		}
		//If there's no process running, wait for the context switch time and set the start and end burst times
		if (!running && readyQueue.size() > 0 && (*readyQueue.begin())->burstStartTime < time) {
			runningProcess = (*readyQueue.begin());
			running = true;
			readyQueue.pop_front();
			runningProcess->burstStartTime = time + csTime / 2;
			runningProcess->burstEndTime = runningProcess->burstStartTime + runningProcess->cpuBurstTimes[runningProcess->burstsCompleted];
		}
		//Adding to each process's wait time and turnaround time when in the ready queue
		for (std::list<process_Info*>::iterator itr = readyQueue.begin(); itr != readyQueue.end(); itr++) {
			(*itr)->waitTimes[(*itr)->burstsCompleted]++;
			(*itr)->turnaroundTimes[(*itr)->burstsCompleted]++;
		}
		time++;
	}



	//Calculate all the statistics
	float avgWaitTime = 0;
	float avgturnAround = 0;
	float avgBurstTime = 0;
	int totalBursts = 0;
	for (int i = 0; i < int(processes.size()); i++) {
		for (int j = 0; j < processes[i].numBursts; j++) {
			avgturnAround += processes[i].turnaroundTimes[j];
			avgWaitTime += processes[i].waitTimes[j];
			avgBurstTime += processes[i].cpuBurstTimes[j];
			totalBursts++;
		}
	}
	avgturnAround /= totalBursts;
	avgWaitTime /= totalBursts;
	avgBurstTime /= totalBursts;
	std::ofstream simout;
	simout.open("simout.txt", std::ios_base::app);
	simout << "Algorithm SJF\n";
	simout << "-- average CPU burst time: " << std::fixed << std::setprecision(3) << avgBurstTime << " ms\n";
	simout << "-- average wait time: " << std::fixed << std::setprecision(3) << avgWaitTime << " ms\n";
	simout << "-- average turnaround time: " << std::fixed << std::setprecision(3) << avgturnAround << " ms\n";
	simout << "-- total number of context switches: " << numContextSwitches << std::endl;
	simout << "-- total number of preemptions: " << numPreemptions << std::endl;
	simout << "-- CPU utilization: " << std::fixed << std::setprecision(3) << cpuUtilization / time * 100 << "%\n";
	simout.close();
}
