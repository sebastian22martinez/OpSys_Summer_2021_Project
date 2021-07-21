#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h> 
#include <sys/shm.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <pthread.h>
#include <math.h>
#include <vector>
#include <iostream>

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
	int burstEndTime = -1;
};

double next_exp(int max, double lambda);
void FCFS(std::vector<process_Info> processes, int csTime);
void SJF(std::vector<process_Info> processes, int csTime);
void SRT(std::vector<process_Info> processes, int csTime);
void RR(std::vector<process_Info> processes, int csTime);

int main(int argc, char ** argv){
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
	srand48(seed);
	std::vector<process_Info> processes;
	for (int i = 0; i < numProc; i++) {
		process_Info process;
		process.processID = 65 + i;
		process.arrivalTime = floor(next_exp(upBound, lambda));
		process.numBursts = ceil(drand48() * 100);
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
}

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

std::string ready_Queue_Format(std::vector<process_Info*> readyQueue) {
	std::string queue = "[Q";
	if (readyQueue.size() == 0) {
		queue += " empty";
	}
	else {
		for (process_Info * p : readyQueue) {
			queue += " ";
			queue += p->processID;
		}
	}
	queue += "]";
	return queue;
}

void FCFS(std::vector<process_Info> processes, int csTime) {
	std::cout << "time 0ms: Simulator started for FCFS [Q empty]" << std::endl;
	int time = 0;
	std::vector<process_Info*> readyQueue;
	process_Info * runningProcess;
	bool running = false;
	int cpuEmptyTime = 0;
	while (true) {
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
		if (readyQueue.size() > 0 && readyQueue[0]->burstStartTime == time) {
			runningProcess = readyQueue[0];
			running = true;
			readyQueue.erase(readyQueue.begin());
			std::cout << "time " << time << "ms: Process " << runningProcess->processID << " started using the CPU for " 
					  << runningProcess->cpuBurstTimes[runningProcess->burstsCompleted] << "ms burst " << ready_Queue_Format(readyQueue) << std::endl;
		}
		if (running && runningProcess->burstEndTime == time) {
			cpuEmptyTime = time + csTime / 2;
			runningProcess->burstsCompleted++;
			std::string burst = "";
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
		for (int i = 0; i < int(processes.size()); i++) {
			if (processes[i].burstEndTime == time) {
				readyQueue.push_back(&processes[i]);
				std::cout << "time " << time << "ms: Process " << processes[i].processID << " completed I/O; added to ready queue " << ready_Queue_Format(readyQueue) << std::endl;
			}
		}
		for (int i = 0; i < int(processes.size()); i++) {
			if (processes[i].arrivalTime == time) {
				readyQueue.push_back(&processes[i]);
				std::cout << "time " << time << "ms: Process " << processes[i].processID << " arrived; added to ready queue " << ready_Queue_Format(readyQueue) << std::endl;
			}
		}
		if (!running && readyQueue.size() > 0 && readyQueue[0]->burstStartTime < time) {
			readyQueue[0]->burstStartTime = time + csTime / 2;
			readyQueue[0]->burstEndTime = readyQueue[0]->burstStartTime + readyQueue[0]->cpuBurstTimes[readyQueue[0]->burstsCompleted];
		}
		time++;
	}
}
