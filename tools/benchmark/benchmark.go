// benchmark.go
// A reusable benchmarking module for Lab Buddy
// Measures execution time and memory usage for any wrapped function

package benchmark

import (
	"fmt"
	"os"
	"runtime"
	"time"
)

// Run wraps any function to measure its runtime and memory usage.
// Additionally reports on host and OS information for repeatability.
func Run(label string, f func()) {
	fmt.Printf("[Benchmark] Running: %s\n", label)

	// Snapshot environment info
	fmt.Println("[Benchmark] Timestamp:", time.Now().Format(time.RFC1123))		// Run begin time
	host, err := os.Hostname()													// Identify hostname 
	if err == nil {
		fmt.Println("[Benchmark] Hostname:", host)								// Report host name
	}
	fmt.Println("[Benchmark] Go Version:", runtime.Version())					// GoLang version
	fmt.Printf("[Benchmark] OS/Arch: %s/%s\n", runtime.GOOS, runtime.GOARCH)	// Operating system

	// Prepare for benchmark
	runtime.GC()																// Measures garbage collection (GC) activity
	var memStart, memEnd runtime.MemStats										// Structs to hold memory statistics before and after execution
	runtime.ReadMemStats(&memStart)												// Capture memory usage before running the function
	start := time.Now()															// Begins running timer
	numCPU := runtime.NumCPU()													// Measures number of available CPUs
	startGoroutines := runtime.NumGoroutine()									// Measures individual Go routines at the beginning of benchmarking

	// Run benchmarked function
	f()																			// Execute the function being benchmarked

	elapsed := time.Since(start)												// Stops running timer
	runtime.ReadMemStats(&memEnd)												// Capture memory usage after the function finishes
	endGoroutines := runtime.NumGoroutine()										// Measures individual Go routines at the end of benchmarking

	// Report resource usage
	fmt.Printf("[Benchmark] Time Elapsed: %v\n", elapsed)																// Reports running time
	fmt.Printf("[Benchmark] Memory Used: %.2f MB\n", float64(memEnd.Alloc-memStart.Alloc)/1024.0/1024.0)				// Shows difference in current heap usage
	fmt.Printf("[Benchmark] Total Allocated: %.2f MB\n", float64(memEnd.TotalAlloc-memStart.TotalAlloc)/1024.0/1024.0)	// Total memory ever allocated during run
	fmt.Printf("[Benchmark] Peak Heap: %.2f MB\n", float64(memEnd.HeapAlloc)/1024.0/1024.0)								// Memory allocated and still in use on the heap after function execution
	fmt.Printf("[Benchmark] GC Cycles: %d\n", memEnd.NumGC-memStart.NumGC)												// Reports GC activity (lower is better)
	fmt.Printf("[Benchmark] Total System Memory Allocated: %.2f MB\n", float64(memEnd.Sys)/1024.0/1024.0)				// Reports all memory requested by the program
	fmt.Printf("[Benchmark] CPU Cores: %d\n", numCPU)																	// Number of available CPU cores
	fmt.Printf("[Benchmark] Goroutines Started: %d → %d\n", startGoroutines, endGoroutines)								// Number of individual Go routines started/ended
	fmt.Println("[Benchmark] ----------------------------------------")													// End
}
