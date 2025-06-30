package sanity_check

import (
	"fmt"
	"lab_buddy_go/config"		// Version control file
)

// Run performs a simple sanity check to ensure Lab_Buddy is
// running properly printing helpful message and version number.
func Run(args []string) {
	fmt.Printf("Successfully running Lab Buddy! (%s)\n", version_control.Main_version)
}